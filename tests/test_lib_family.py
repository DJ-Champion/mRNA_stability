"""Unit tests for lib/family.py.

Emphasis is on the three properties whose failure would be silent and would
corrupt the CV blocking rather than crash a run:

  1. Translation must not truncate at internal stops (the seqkit --trim trap).
  2. Singletons must survive clustering.
  3. Family labels must be reproducible across runs.
"""
import pytest

from lib.family import (
    Level, SearchParams, UnionFind, canonical_families, iter_fasta,
    make_protein_id, normalised_bitscore, post_process_translation,
    search_hash, split_protein_id, validate_level_within_search,
)


# ---------------------------------------------------------------------------
# Translation post-processing
# ---------------------------------------------------------------------------

class TestPostProcessTranslation:
    """Regression guard against `seqkit translate --trim`.

    The expected values are the correct translations; --trim (verified
    against seqkit v2.13.0) instead returns 'MAWK', 'MA', 'MA' for these
    three inputs, silently reducing the last two to N-terminal stubs.
    """

    def test_terminal_stop_stripped(self):
        res = post_process_translation('MAWK*')
        assert res.aa == 'MAWK'
        assert res.had_terminal_stop is True
        assert res.had_internal_stop is False

    def test_internal_stop_becomes_X_not_truncation(self):
        res = post_process_translation('MA*KA*')
        assert res.aa == 'MAXKA'          # --trim would give 'MA'
        assert res.had_internal_stop is True
        assert res.had_terminal_stop is True

    def test_ambiguous_residue_preserved(self):
        res = post_process_translation('MAXK*')
        assert res.aa == 'MAXK'           # --trim would give 'MA'
        assert res.had_internal_stop is False

    def test_only_one_terminal_stop_stripped(self):
        # Two trailing stops: the second is internal once the first is gone.
        res = post_process_translation('MAWK**')
        assert res.aa == 'MAWKX'
        assert res.had_internal_stop is True

    def test_no_terminal_stop_flagged(self):
        res = post_process_translation('MAWK')
        assert res.aa == 'MAWK'
        assert res.had_terminal_stop is False

    def test_selenoprotein_keeps_full_length(self):
        # A UGA-readthrough protein must not lose its C-terminal domain.
        aa = 'M' + 'A' * 100 + '*' + 'K' * 100 + '*'
        res = post_process_translation(aa)
        assert len(res.aa) == 202
        assert res.aa.count('X') == 1


# ---------------------------------------------------------------------------
# Normalised bitscore
# ---------------------------------------------------------------------------

class TestNormalisedBitscore:

    def test_self_comparison_is_one(self):
        assert normalised_bitscore(500.0, 500.0, 500.0) == pytest.approx(1.0)

    def test_uses_max_denominator_to_reject_domain_match(self):
        # 100 aa protein fully matching one domain of a 900 aa protein.
        nbs = normalised_bitscore(bits=200.0, self_bits_q=200.0, self_bits_t=1800.0)
        assert nbs == pytest.approx(200 / 1800, abs=1e-6)
        assert nbs < 0.25          # rejected even by the loosest level

    def test_min_denominator_would_have_chained(self):
        # Documents why max() is used: min() accepts the same edge.
        assert 200.0 / min(200.0, 1800.0) == pytest.approx(1.0)

    def test_symmetric(self):
        a = normalised_bitscore(150.0, 200.0, 1800.0)
        b = normalised_bitscore(150.0, 1800.0, 200.0)
        assert a == b

    def test_clamped_to_one(self):
        assert normalised_bitscore(600.0, 500.0, 400.0) == 1.0

    def test_zero_self_bits_is_zero_not_error(self):
        assert normalised_bitscore(10.0, 0.0, 0.0) == 0.0


# ---------------------------------------------------------------------------
# Level filtering
# ---------------------------------------------------------------------------

class TestLevel:

    @staticmethod
    def _level():
        return Level('medium', min_nbs=0.40, min_cov=0.65,
                     max_evalue=1e-5, min_fident=0.30)

    def test_accepts_when_all_criteria_met(self):
        assert self._level().accepts(0.5, 0.7, 0.7, 1e-8, 0.4)

    def test_requires_both_coverages(self):
        lv = self._level()
        # Query well covered, target barely — the domain-in-larger-protein
        # geometry. cov_mode 0 semantics mean this must be rejected.
        assert not lv.accepts(0.5, 0.95, 0.20, 1e-8, 0.4)
        assert not lv.accepts(0.5, 0.20, 0.95, 1e-8, 0.4)

    @pytest.mark.parametrize('nbs,qcov,tcov,ev,fid', [
        (0.39, 0.70, 0.70, 1e-8, 0.40),   # nbs below
        (0.50, 0.64, 0.70, 1e-8, 0.40),   # qcov below
        (0.50, 0.70, 0.70, 1e-4, 0.40),   # evalue above
        (0.50, 0.70, 0.70, 1e-8, 0.29),   # fident below
    ])
    def test_rejects_on_any_single_criterion(self, nbs, qcov, tcov, ev, fid):
        assert not self._level().accepts(nbs, qcov, tcov, ev, fid)


class TestValidateLevelWithinSearch:

    @staticmethod
    def _search():
        return SearchParams(sensitivity=7.5, max_seqs=2000, evalue=1e-3,
                            min_cov=0.4, cov_mode=0, threads=4)

    def test_tighter_level_passes(self):
        lv = Level('strict', 0.6, 0.8, 1e-10, 0.5)
        assert validate_level_within_search(lv, self._search()) == []

    def test_equal_to_search_passes(self):
        lv = Level('loose', 0.25, 0.4, 1e-3, 0.2)
        assert validate_level_within_search(lv, self._search()) == []

    def test_looser_evalue_rejected(self):
        lv = Level('bad', 0.25, 0.5, 1e-1, 0.2)
        assert any('evalue' in p for p in validate_level_within_search(lv, self._search()))

    def test_looser_coverage_rejected(self):
        lv = Level('bad', 0.25, 0.2, 1e-5, 0.2)
        assert any('min_cov' in p for p in validate_level_within_search(lv, self._search()))


# ---------------------------------------------------------------------------
# Union-find
# ---------------------------------------------------------------------------

class TestUnionFind:

    def test_starts_fully_disjoint(self):
        uf = UnionFind(5)
        assert len({uf.find(i) for i in range(5)}) == 5

    def test_transitive_closure(self):
        # Single linkage: a~b and b~c puts a and c together even though the
        # a-c edge was never observed. This is the leakage guarantee.
        uf = UnionFind(3)
        uf.union(0, 1)
        uf.union(1, 2)
        assert uf.find(0) == uf.find(2)

    def test_union_reports_whether_merge_happened(self):
        uf = UnionFind(3)
        assert uf.union(0, 1) is True
        assert uf.union(1, 0) is False

    def test_long_chain(self):
        uf = UnionFind(1000)
        for i in range(999):
            uf.union(i, i + 1)
        assert len({uf.find(i) for i in range(1000)}) == 1


# ---------------------------------------------------------------------------
# Canonical family assignment
# ---------------------------------------------------------------------------

class TestCanonicalFamilies:

    def test_singletons_are_included(self):
        # The load-bearing property: genes with no edges must still get a
        # family. Building the graph from the edge list alone drops them.
        ids = ['a', 'b', 'c', 'd']
        uf = UnionFind(4)
        uf.union(0, 1)
        labels, groups = canonical_families(uf, ids, 'medium')
        assert all(lbl for lbl in labels)
        assert len(labels) == 4
        assert len(groups) == 3          # {a,b}, {c}, {d}

    def test_ordered_by_size_descending(self):
        ids = list('abcdef')
        uf = UnionFind(6)
        uf.union(0, 1); uf.union(1, 2)   # 3-member family
        uf.union(3, 4)                   # 2-member family
        labels, groups = canonical_families(uf, ids, 'medium')
        assert [len(g) for g in groups] == [3, 2, 1]
        assert labels[0] == 'fam_medium_00001'
        assert labels[3] == 'fam_medium_00002'
        assert labels[5] == 'fam_medium_00003'

    def test_ties_broken_by_smallest_member_id(self):
        ids = ['zebra', 'yak', 'ant', 'bee']
        uf = UnionFind(4)
        uf.union(0, 1)                   # {zebra, yak}
        uf.union(2, 3)                   # {ant, bee}  -> smaller min id
        labels, _ = canonical_families(uf, ids, 'x')
        assert labels[2] == 'fam_x_00001'
        assert labels[0] == 'fam_x_00002'

    def test_labels_independent_of_union_order(self):
        # mmseqs is threaded; edge arrival order must not change the labels.
        ids = list('abcdef')
        uf1 = UnionFind(6)
        for a, b in [(0, 1), (1, 2), (3, 4)]:
            uf1.union(a, b)
        uf2 = UnionFind(6)
        for a, b in [(4, 3), (2, 1), (1, 0)]:
            uf2.union(a, b)
        assert canonical_families(uf1, ids, 'm')[0] == \
               canonical_families(uf2, ids, 'm')[0]

    def test_every_index_labelled(self):
        ids = [f's{i}' for i in range(50)]
        uf = UnionFind(50)
        uf.union(0, 49)
        labels, groups = canonical_families(uf, ids, 'loose')
        assert sum(len(g) for g in groups) == 50
        assert '' not in labels


# ---------------------------------------------------------------------------
# Protein IDs
# ---------------------------------------------------------------------------

class TestProteinIds:

    def test_round_trip(self):
        pid = make_protein_id('human', 'ENSG00000141510', 'ENST00000269305.9')
        assert pid == 'human|ENSG00000141510|ENST00000269305.9'
        assert split_protein_id(pid) == (
            'human', 'ENSG00000141510', 'ENST00000269305.9')

    @pytest.mark.parametrize('species,gene,tx', [
        ('hu|man', 'G1', 'T1'),
        ('human', 'G|1', 'T1'),
        ('human', 'G1', 'T 1'),
        ('home sapiens', 'G1', 'T1'),
    ])
    def test_rejects_ambiguous_components(self, species, gene, tx):
        with pytest.raises(ValueError):
            make_protein_id(species, gene, tx)

    def test_split_rejects_malformed(self):
        with pytest.raises(ValueError):
            split_protein_id('human|ENSG0001')


# ---------------------------------------------------------------------------
# FASTA reading
# ---------------------------------------------------------------------------

class TestIterFasta:

    def test_reads_wrapped_sequences(self, tmp_path):
        # seqkit wraps at 60 by default.
        p = tmp_path / 'w.faa'
        p.write_text('>one\nMAWK\n>two\n' + 'A' * 60 + '\n' + 'K' * 10 + '\n')
        recs = dict(iter_fasta(p))
        assert recs['one'] == 'MAWK'
        assert recs['two'] == 'A' * 60 + 'K' * 10

    def test_header_truncated_at_whitespace(self, tmp_path):
        p = tmp_path / 'd.faa'
        p.write_text('>id1 some description here\nMAWK\n')
        assert list(iter_fasta(p)) == [('id1', 'MAWK')]

    def test_empty_file(self, tmp_path):
        p = tmp_path / 'e.faa'
        p.write_text('')
        assert list(iter_fasta(p)) == []


# ---------------------------------------------------------------------------
# Search hash
# ---------------------------------------------------------------------------

class TestSearchHash:

    @staticmethod
    def _search(**kw):
        base = dict(sensitivity=7.5, max_seqs=2000, evalue=1e-3,
                    min_cov=0.4, cov_mode=0, threads=12)
        base.update(kw)
        return SearchParams(**base)

    def test_stable_across_calls(self):
        m = [('human_test', 'human', 1)]
        assert search_hash(self._search(), m) == search_hash(self._search(), m)

    def test_member_order_irrelevant(self):
        a = [('a', 'human', 1), ('b', 'mouse', 1)]
        assert search_hash(self._search(), a) == \
               search_hash(self._search(), list(reversed(a)))

    def test_threads_do_not_change_identity(self):
        # Thread count must not change the identity of the result.
        m = [('human_test', 'human', 1)]
        assert search_hash(self._search(threads=1), m) == \
               search_hash(self._search(threads=64), m)

    def test_search_params_change_hash(self):
        m = [('human_test', 'human', 1)]
        assert search_hash(self._search(), m) != \
               search_hash(self._search(evalue=1e-5), m)

    def test_added_member_changes_hash(self):
        assert search_hash(self._search(), [('a', 'human', 1)]) != \
               search_hash(self._search(), [('a', 'human', 1), ('b', 'mouse', 1)])
