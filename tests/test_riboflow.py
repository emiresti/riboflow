"""Comprehensive test suite for RiboFlow."""

from __future__ import annotations

import unittest

from riboflow import (
    ORGANISMS,
    SAMPLES,
    analyze_sequence,
    apply_synonymous_mutations,
    compare_harmonization,
    get_organism,
    harmonize_codons,
    naive_optimize_codons,
    parse_fasta_or_raw,
    scan_critical_silent_mutations,
    simulate_ribosome_traffic,
    translate_codons,
    validate_cds,
)
from riboflow.web.server import get_metadata_payload, handle_api_post


class TestOrganismsAndSamples(unittest.TestCase):
    def test_all_organisms_have_64_codons(self) -> None:
        for key, org in ORGANISMS.items():
            self.assertEqual(len(org.codon_freq_per_1000), 64, msg=key)
            rates = org.elongation_rates()
            self.assertEqual(len(rates), 64)
            self.assertTrue(all(r > 0 for r in rates.values()))

    def test_all_curated_samples_are_valid_cds(self) -> None:
        for sample_id, sample in SAMPLES.items():
            codons = validate_cds(sample.sequence)
            self.assertGreaterEqual(len(codons), 50, msg=sample_id)
            self.assertEqual(codons[0], "ATG", msg=sample_id)
            self.assertIn(codons[-1], {"TAA", "TAG", "TGA"}, msg=sample_id)


class TestSequenceAnalyzer(unittest.TestCase):
    def test_fasta_and_rna_normalization(self) -> None:
        raw = ">TestGene\nAUGGCUGCUAAAGCUGCUGCUAAAGCUGCUGCUAAAGCUGCUGCUAAAUAA\n"
        header, dna = parse_fasta_or_raw(raw)
        self.assertEqual(header, "TestGene")
        self.assertNotIn("U", dna)
        codons = validate_cds(dna)
        self.assertEqual(codons[0], "ATG")

    def test_invalid_cds_rejected(self) -> None:
        with self.assertRaises(ValueError):
            validate_cds("ATGCG")  # not multiple of 3
        with self.assertRaises(ValueError):
            validate_cds("GCT" * 20)  # missing ATG start
        with self.assertRaises(ValueError):
            validate_cds("ATG" + "GCT" * 10 + "TAA" + "GCT" * 10)  # internal stop

    def test_analysis_metrics_bounds(self) -> None:
        ecoli = get_organism("ecoli")
        res = analyze_sequence(SAMPLES["human_ins"].sequence, ecoli)
        self.assertGreater(res.cai, 0.0)
        self.assertLessEqual(res.cai, 1.0)
        self.assertGreater(res.gc_percent, 40.0)
        self.assertEqual(len(res.codons), res.codon_count)


class TestTASEPSimulator(unittest.TestCase):
    def test_steric_exclusion_never_violated(self) -> None:
        ecoli = get_organism("ecoli")
        sim = simulate_ribosome_traffic(
            SAMPLES["jam_benchmark"].sequence,
            ecoli,
            initiation_rate=1.6,
            duration_seconds=80.0,
            capture_frames=80,
            seed=7,
        )
        self.assertGreater(sim.completed_proteins, 0)
        self.assertGreater(sim.collision_fraction, 0.0)
        for frame in sim.frames:
            positions = sorted(r.pos for r in frame.ribosomes)
            for a, b in zip(positions, positions[1:]):
                self.assertGreaterEqual(b - a, ecoli.ribosome_footprint)

    def test_deterministic_with_seed(self) -> None:
        human = get_organism("human")
        s1 = simulate_ribosome_traffic(SAMPLES["cftr_nbd1"].sequence, human, seed=123)
        s2 = simulate_ribosome_traffic(SAMPLES["cftr_nbd1"].sequence, human, seed=123)
        self.assertEqual(s1.completed_proteins, s2.completed_proteins)
        self.assertEqual(s1.collision_events, s2.collision_events)
        self.assertEqual(s1.codon_occupancy, s2.codon_occupancy)

    def test_bottleneck_causes_more_collisions_than_fast_optimized(self) -> None:
        ecoli = get_organism("ecoli")
        jam_seq = SAMPLES["jam_benchmark"].sequence
        fast_seq = "".join(naive_optimize_codons(validate_cds(jam_seq), ecoli))

        jam_sim = simulate_ribosome_traffic(jam_seq, ecoli, initiation_rate=1.4, seed=42)
        fast_sim = simulate_ribosome_traffic(fast_seq, ecoli, initiation_rate=1.4, seed=42)

        self.assertGreater(jam_sim.ngd_risk_score, fast_sim.ngd_risk_score)
        self.assertGreater(jam_sim.collision_fraction, fast_sim.collision_fraction)


class TestCodonHarmonizer(unittest.TestCase):
    def test_harmonization_preserves_exact_protein_and_improves_rhythm(self) -> None:
        human = get_organism("human")
        ecoli = get_organism("ecoli")
        seq = SAMPLES["human_ins"].sequence

        report = compare_harmonization(seq, source=human, host=ecoli, duration_seconds=60.0)
        self.assertTrue(report.protein_preserved)
        self.assertEqual(
            translate_codons(validate_cds(report.harmonized.sequence)),
            report.protein_sequence,
        )
        # Harmonized rhythm should correlate better with native source rhythm than naive CAI max
        self.assertGreater(
            report.harmonized.rhythm_correlation,
            report.naive_optimized.rhythm_correlation,
        )

    def test_harmonize_codons_handles_all_samples(self) -> None:
        for sample in SAMPLES.values():
            src = get_organism(sample.organism)
            host = get_organism(sample.recommended_host)
            wt = validate_cds(sample.sequence)
            harm = harmonize_codons(wt, src, host)
            self.assertEqual(translate_codons(wt), translate_codons(harm))


class TestSilentMutationLab(unittest.TestCase):
    def test_synonymous_mutation_enforced(self) -> None:
        ecoli = get_organism("ecoli")
        seq = SAMPLES["human_ins"].sequence
        codons = validate_cds(seq)
        # Find a Leucine or Arginine codon to mutate silently
        target_idx = next(i for i in range(5, len(codons) - 1) if codons[i] == "CTG")
        comp = apply_synonymous_mutations(seq, ecoli, edits={target_idx: "CTA"}, duration_seconds=60.0)
        self.assertTrue(comp.protein_identical)
        self.assertEqual(len(comp.edits), 1)

        # Non-synonymous edit must raise ValueError
        with self.assertRaises(ValueError):
            apply_synonymous_mutations(seq, ecoli, edits={target_idx: "AAA"})

    def test_hotspot_scanner_returns_valid_silent_swaps(self) -> None:
        ecoli = get_organism("ecoli")
        hotspots = scan_critical_silent_mutations(
            SAMPLES["cftr_nbd1"].sequence,
            ecoli,
            initiation_rate=1.3,
            top_k=4,
            seed=42,
        )
        self.assertGreater(len(hotspots), 0)
        for h in hotspots:
            self.assertEqual(h.nucleotide_changes, 1)
            self.assertGreater(h.rate_drop_percent, 0.0)


class TestWebAPIRouter(unittest.TestCase):
    def test_meta_and_api_endpoints(self) -> None:
        meta = get_metadata_payload()
        self.assertEqual(len(meta["organisms"]), 4)
        self.assertEqual(len(meta["samples"]), 4)

        status, body = handle_api_post(
            "/api/simulate",
            {"sequence": SAMPLES["human_ins"].sequence, "organism": "ecoli", "duration": 40},
        )
        self.assertEqual(status, 200)
        self.assertIn("analysis", body)
        self.assertIn("simulation", body)

        status_h, body_h = handle_api_post(
            "/api/harmonize",
            {"sequence": SAMPLES["human_ins"].sequence, "source": "human", "host": "ecoli"},
        )
        self.assertEqual(status_h, 200)
        self.assertTrue(body_h["protein_preserved"])

        status_err, body_err = handle_api_post(
            "/api/simulate",
            {"sequence": "NOT_DNA_123", "organism": "ecoli"},
        )
        self.assertEqual(status_err, 400)
        self.assertIn("error", body_err)


class TestCLI(unittest.TestCase):
    def test_cli_subcommands_run_cleanly(self) -> None:
        from riboflow.cli import main

        self.assertEqual(main(["analyze", "--sample", "human_ins", "-o", "ecoli"]), 0)
        self.assertEqual(main(["simulate", "--sample", "jam_benchmark", "-o", "ecoli", "-d", "40"]), 0)
        self.assertEqual(main(["harmonize", "--sample", "human_ins", "-s", "human", "-t", "ecoli", "-d", "40"]), 0)
        self.assertEqual(main(["mutate", "--sample", "cftr_nbd1", "-o", "ecoli", "--scan", "--top-k", "3"]), 0)


if __name__ == "__main__":
    unittest.main()
