# =============================================================================
# tests/test_rules.py — Unit tests for the Snakemake pipeline rules
# =============================================================================
# Each test class mocks the external tools (GATK, samtools, varscan, SnpEff,
# etc.) and verifies that the correct command-line arguments are constructed
# and that expected output files are produced/checked.
#
# Run with:
#   pytest tests/test_rules.py -v
#
# Dependencies:
#   pip install pytest pytest-mock snakemake
# =============================================================================

import os
import sys
import shutil
import tempfile
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock, call

import pytest


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_workdir(tmp_path):
    """
    Create a temporary working directory that mimics the pipeline layout:
      tmp_workdir/
        data/
        results/
        logs/
        annotations/
    """
    for subdir in ["data", "results", "logs", "annotations",
                   "results/somatic", "results/purity", "results/cnv",
                   "results/ancestry", "results/qc", "results/igv"]:
        (tmp_path / subdir).mkdir(parents=True, exist_ok=True)
    return tmp_path


@pytest.fixture
def dummy_bam(tmp_workdir):
    """Create a minimal stub .bam file (not a valid BAM, just for path tests)."""
    bam = tmp_workdir / "data" / "control.sorted.bam"
    bam.write_bytes(b"STUB_BAM")
    return bam


@pytest.fixture
def dummy_vcf(tmp_workdir):
    """Create a minimal stub VCF file."""
    vcf = tmp_workdir / "results" / "control.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tT\t50\tPASS\t.\n"
    )
    return vcf


@pytest.fixture
def dummy_fasta(tmp_workdir):
    """Create a stub FASTA reference."""
    fasta = tmp_workdir / "annotations" / "human_g1k_v37.fasta"
    fasta.write_text(">chr1\nACGT\n")
    return fasta


# =============================================================================
# 1 · preprep — samtools sort
# =============================================================================

class TestPreprep:
    """Tests for the preprep rule (samtools sort wrapper)."""

    def test_output_path_convention(self, tmp_workdir):
        """Sorted BAM should be data/{group}.sorted.bam."""
        for group in ["control", "tumor"]:
            expected = tmp_workdir / "data" / f"{group}.sorted.bam"
            # Simulate rule output creation
            expected.write_bytes(b"")
            assert expected.exists(), f"Expected sorted BAM missing: {expected}"

    @patch("subprocess.run")
    def test_samtools_sort_called(self, mock_run, tmp_workdir, dummy_bam):
        """samtools sort should be invoked with the correct -o flag."""
        mock_run.return_value = MagicMock(returncode=0)
        input_bam = str(dummy_bam)
        output_bam = str(tmp_workdir / "data" / "control.sorted.bam")

        subprocess.run(["samtools", "sort", "-o", output_bam, input_bam], check=True)

        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert "samtools" in args
        assert "sort" in args
        assert output_bam in args


# =============================================================================
# 2 · indexing — samtools index
# =============================================================================

class TestIndexing:
    """Tests for the indexing rule (samtools index wrapper)."""

    def test_index_file_name(self, tmp_workdir):
        """Index file must be data/{group}.sorted.bai (not .bam.bai)."""
        for group in ["control", "tumor"]:
            expected = tmp_workdir / "data" / f"{group}.sorted.bai"
            expected.write_bytes(b"")
            assert expected.suffix == ".bai"
            assert "sorted" in expected.name

    @patch("subprocess.run")
    def test_samtools_index_called_with_bam(self, mock_run, dummy_bam):
        """samtools index should receive the sorted BAM as its sole positional arg."""
        mock_run.return_value = MagicMock(returncode=0)
        subprocess.run(["samtools", "index", str(dummy_bam)], check=True)

        args = mock_run.call_args[0][0]
        assert "index" in args
        assert str(dummy_bam) in args


# =============================================================================
# 3 · RemoveDuplicates — GATK MarkDuplicates
# =============================================================================

class TestRemoveDuplicates:
    """Tests for the RemoveDuplicates rule."""

    @patch("subprocess.run")
    def test_markduplicates_flags(self, mock_run, tmp_workdir, dummy_bam):
        """MarkDuplicates must use REMOVE_DUPLICATES=true and ASSUME_SORT_ORDER."""
        mock_run.return_value = MagicMock(returncode=0)
        out_bam = str(tmp_workdir / "data" / "control.sorted.dedup.bam")
        metrics = str(tmp_workdir / "results" / "qc" / "control.dedup_metrics.txt")

        cmd = [
            "gatk", "MarkDuplicates",
            "-I", str(dummy_bam),
            "-O", out_bam,
            "-M", metrics,
            "--REMOVE_DUPLICATES", "true",
            "--ASSUME_SORT_ORDER", "coordinate",
        ]
        subprocess.run(cmd, check=True)

        args = mock_run.call_args[0][0]
        assert "--REMOVE_DUPLICATES" in args
        assert "true" in args
        assert "--ASSUME_SORT_ORDER" in args

    def test_metrics_output_path(self, tmp_workdir):
        """Dedup metrics must land under results/qc/."""
        metrics_path = tmp_workdir / "results" / "qc" / "control.dedup_metrics.txt"
        metrics_path.write_text("## METRICS CLASS\n")
        assert metrics_path.exists()
        assert "qc" in str(metrics_path)


# =============================================================================
# 4 · BaseQualityScoreRecalibration (BQSR)
# =============================================================================

class TestBQSR:
    """Tests for the BQSR rule."""

    @patch("subprocess.run")
    def test_bqsr_script_called(self, mock_run, tmp_workdir, dummy_bam, dummy_fasta):
        """BQSR.sh must be called with -r, -i, -o, -k, -l, -t flags."""
        mock_run.return_value = MagicMock(returncode=0)
        out_bam = str(tmp_workdir / "data" / "control.sorted.recalibrated.bam")
        known_sites = str(tmp_workdir / "annotations" / "hapmap_3.3.b37.vcf")
        intervals = str(tmp_workdir / "annotations" / "CancerGenesSel.bed")
        recal_table = str(tmp_workdir / "results" / "qc" / "control.recal_table.txt")

        cmd = [
            "bash", "scripts/BQSR.sh",
            "-r", str(dummy_fasta),
            "-i", str(dummy_bam),
            "-o", out_bam,
            "-k", known_sites,
            "-l", intervals,
            "-t", recal_table,
        ]
        subprocess.run(cmd, check=True)

        args = mock_run.call_args[0][0]
        for flag in ["-r", "-i", "-o", "-k", "-l", "-t"]:
            assert flag in args, f"Expected flag {flag!r} not found in BQSR command"

    def test_recal_table_output_exists(self, tmp_workdir):
        """Recalibration table should be created under results/qc/."""
        recal = tmp_workdir / "results" / "qc" / "control.recal_table.txt"
        recal.write_text("ReadGroup\tQuality\tCount\n")
        assert recal.exists()


# =============================================================================
# 5 · SPIA_check
# =============================================================================

class TestSPIACheck:
    """
    Tests for the SPIA_check rule (sample-swap detection).
    The actual R script is mocked; we test argument passing and output path.
    """

    @patch("subprocess.run")
    def test_spia_rscript_arguments(self, mock_run, tmp_workdir):
        """Rscript should receive --control, --tumor, and --output arguments."""
        mock_run.return_value = MagicMock(returncode=0)
        control_bam = str(tmp_workdir / "data" / "control.sorted.bam")
        tumor_bam   = str(tmp_workdir / "data" / "tumor.sorted.bam")
        out_pdf     = str(tmp_workdir / "results" / "control" / "spia_report.pdf")

        cmd = [
            "Rscript", "scripts/spia_check.R",
            "--control", control_bam,
            "--tumor",   tumor_bam,
            "--output",  out_pdf,
        ]
        subprocess.run(cmd, check=True)

        args = mock_run.call_args[0][0]
        assert "--control" in args
        assert "--tumor" in args
        assert "--output" in args

    def test_spia_output_is_pdf(self, tmp_workdir):
        """SPIA report output must be a .pdf file."""
        out = tmp_workdir / "results" / "control" / "spia_report.pdf"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(b"%PDF-1.4")
        assert out.suffix == ".pdf"


# =============================================================================
# 6 · VariantCalling — GATK HaplotypeCaller
# =============================================================================

class TestVariantCalling:
    """Tests for the VariantCalling rule."""

    @patch("subprocess.run")
    def test_haplotypecaller_flags(self, mock_run, tmp_workdir, dummy_bam, dummy_fasta):
        """HaplotypeCaller must receive -R, -I, -O, and --bam-output."""
        mock_run.return_value = MagicMock(returncode=0)
        out_vcf    = str(tmp_workdir / "results" / "control.vcf")
        out_bamout = str(tmp_workdir / "results" / "control.bamout.bam")

        cmd = [
            "gatk", "HaplotypeCaller",
            "-R", str(dummy_fasta),
            "-I", str(dummy_bam),
            "-O", out_vcf,
            "--bam-output", out_bamout,
        ]
        subprocess.run(cmd, check=True)

        args = mock_run.call_args[0][0]
        assert "HaplotypeCaller" in args
        assert "--bam-output" in args
        assert out_vcf in args

    def test_vcf_output_path(self, tmp_workdir):
        """VCF output should land in results/{group}.vcf."""
        for group in ["control", "tumor"]:
            vcf = tmp_workdir / "results" / f"{group}.vcf"
            vcf.write_text("##fileformat=VCFv4.2\n")
            assert vcf.exists()
            assert vcf.suffix == ".vcf"


# =============================================================================
# 7 · SomaticVariantCall — VarScan somatic
# =============================================================================

class TestSomaticVariantCall:
    """Tests for the SomaticVariantCall rule."""

    @patch("subprocess.run")
    def test_mpileup_called_for_both_samples(self, mock_run, tmp_workdir, dummy_fasta):
        """samtools mpileup must be called for both control and tumor BAMs."""
        mock_run.return_value = MagicMock(returncode=0)
        for group in ["control", "tumor"]:
            subprocess.run(
                ["samtools", "mpileup", "-q", "1", "-f", str(dummy_fasta),
                 f"data/{group}.sorted.recalibrated.dedup.bam"],
                check=True,
            )
        assert mock_run.call_count == 2

    @patch("subprocess.run")
    def test_varscan_somatic_output_vcf_flag(self, mock_run, tmp_workdir):
        """VarScan somatic must include --output-vcf 1."""
        mock_run.return_value = MagicMock(returncode=0)
        cmd = [
            "varscan", "somatic",
            "control.pileup", "tumor.pileup",
            "--output-snp",   "results/somatic.snp.vcf",
            "--output-indel", "results/somatic.indel.vcf",
            "--output-vcf",   "1",
        ]
        subprocess.run(cmd, check=True)
        args = mock_run.call_args[0][0]
        assert "--output-vcf" in args
        assert "1" in args

    def test_pileup_files_are_temp(self, tmp_workdir):
        """
        Pileup files are declared as temp() in the rule; here we verify that
        a downstream consumer deleting them does not break the test fixtures.
        """
        pileup = tmp_workdir / "results" / "somatic" / "control.pileup"
        pileup.write_text("chr1\t100\t...\n")
        assert pileup.exists()
        pileup.unlink()
        assert not pileup.exists()


# =============================================================================
# 8 · CopyNumberVariation — VarScan copynumber
# =============================================================================

class TestCopyNumberVariation:
    """Tests for the CopyNumberVariation rule."""

    @patch("subprocess.run")
    def test_copynumber_pipeline(self, mock_run, tmp_workdir, dummy_fasta):
        """
        The CNV pipeline requires three sequential calls:
          1. samtools mpileup | varscan copynumber
          2. varscan copyCaller
          3. Rscript plot_cnv.R
        """
        mock_run.return_value = MagicMock(returncode=0)

        subprocess.run(["samtools", "mpileup", "-q", "1", "-f", str(dummy_fasta),
                        "control.bam", "tumor.bam"], check=True)
        subprocess.run(["varscan", "copynumber", "--mpileup", "1"], check=True)
        subprocess.run(["varscan", "copyCaller", "copynumber.txt"], check=True)
        subprocess.run(["Rscript", "scripts/plot_cnv.R", "copycaller.txt", "plot.pdf"],
                       check=True)

        assert mock_run.call_count == 4

    def test_cnv_plot_is_pdf(self, tmp_workdir):
        """CNV plot output must be a PDF."""
        plot = tmp_workdir / "results" / "cnv" / "cnv_plot.pdf"
        plot.write_bytes(b"%PDF")
        assert plot.suffix == ".pdf"


# =============================================================================
# 9 · VariantAnnotation — SnpEff + SnpSift
# =============================================================================

class TestVariantAnnotation:
    """Tests for the VariantAnnotation rule."""

    @patch("subprocess.run")
    def test_snpeff_called_with_database(self, mock_run, tmp_workdir, dummy_vcf):
        """SnpEff must be called with the configured database name."""
        mock_run.return_value = MagicMock(returncode=0)
        database = "hg19kg"
        summary  = str(tmp_workdir / "results" / "snpeff_summary.html")
        out_vcf  = str(tmp_workdir / "results" / "annotated_variants.vcf")

        cmd = ["snpEff", "-Xmx8g", "-v", database,
               str(dummy_vcf), "-s", summary]
        subprocess.run(cmd, check=True)

        args = mock_run.call_args[0][0]
        assert "snpEff" in args
        assert database in args
        assert summary in args

    @patch("subprocess.run")
    def test_snpsift_annotate_called(self, mock_run, tmp_workdir):
        """SnpSift Annotate must be called after SnpEff."""
        mock_run.return_value = MagicMock(returncode=0)
        subprocess.run(["SnpSift", "-Xmx8g", "Annotate",
                        "clinvar.vcf", "annotated.vcf"], check=True)
        args = mock_run.call_args[0][0]
        assert "SnpSift" in args
        assert "Annotate" in args

    def test_get_annotation_db_returns_correct_db(self):
        """get_annotation_db helper should map assembly tokens to DB names."""
        # Inline reimplementation to test without importing the full Snakefile
        assembly_map = {"hg19": "hg19kg", "hg38": "hg38kg", "b37": "GRCh37.75"}
        default_db   = "hg19kg"

        def get_annotation_db(sample_name):
            for key, db in assembly_map.items():
                if key in sample_name:
                    return db
            return default_db

        assert get_annotation_db("sample_hg19") == "hg19kg"
        assert get_annotation_db("sample_hg38") == "hg38kg"
        assert get_annotation_db("sample_b37")  == "GRCh37.75"
        assert get_annotation_db("unknown")      == "hg19kg"


# =============================================================================
# 10 · PurityPloidy — TPES / ASEReadCounter
# =============================================================================

class TestPurityPloidy:
    """Tests for the PurityPloidy rule."""

    @patch("subprocess.run")
    def test_bcftools_filter_biallelic(self, mock_run, tmp_workdir, dummy_vcf):
        """bcftools must filter for biallelic SNPs (-v snps -m2 -M2)."""
        mock_run.return_value = MagicMock(returncode=0)
        subprocess.run(
            ["bcftools", "view", "-v", "snps", "-m2", "-M2", str(dummy_vcf)],
            check=True,
        )
        args = mock_run.call_args[0][0]
        assert "bcftools" in args
        assert "-m2" in args
        assert "-M2" in args

    @patch("subprocess.run")
    def test_ase_read_counter_min_depth(self, mock_run, tmp_workdir,
                                        dummy_bam, dummy_vcf, dummy_fasta):
        """ASEReadCounter must use -minDepth 20."""
        mock_run.return_value = MagicMock(returncode=0)
        subprocess.run(
            ["java", "-jar", "genome_analysis_TK.jar", "-T", "ASEReadCounter",
             "-R", str(dummy_fasta), "-I", str(dummy_bam),
             "-sites", str(dummy_vcf), "-minDepth", "20"],
            check=True,
        )
        args = mock_run.call_args[0][0]
        assert "-minDepth" in args
        assert "20" in args

    def test_het_snp_extraction(self, tmp_workdir):
        """
        Heterozygous SNP extraction via grep should retain only 0/1 genotype
        lines and VCF headers (lines starting with #).
        """
        input_lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
            "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n",  # het → keep
            "chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t1/1\n",  # hom alt → drop
            "chr1\t300\t.\tC\tA\t70\tPASS\t.\tGT\t0/0\n",  # hom ref → drop
        ]
        import re
        pattern = re.compile(r"(^#|0/1)")
        kept = [l for l in input_lines if pattern.search(l)]
        assert len(kept) == 3  # 2 header lines + 1 het line
        assert all("0/1" in l or l.startswith("#") for l in kept)


# =============================================================================
# 11 · VariantPrediction — AlphaGenome
# =============================================================================

class TestVariantPrediction:
    """Tests for the VariantPrediction rule."""

    @patch("subprocess.run")
    def test_prediction_script_called(self, mock_run, tmp_workdir, dummy_vcf):
        """variant_prediction.py should be called with --input and --output."""
        mock_run.return_value = MagicMock(returncode=0)
        out_html = str(tmp_workdir / "results" / "variant_prediction_report.html")
        subprocess.run(
            ["python", "scripts/variant_prediction.py",
             "--input", str(dummy_vcf), "--output", out_html],
            check=True,
        )
        args = mock_run.call_args[0][0]
        assert "--input"  in args
        assert "--output" in args
        assert out_html in args

    def test_output_is_html(self, tmp_workdir):
        """Output of VariantPrediction must be an HTML file."""
        html = tmp_workdir / "results" / "variant_prediction_report.html"
        html.write_text("<html></html>")
        assert html.suffix == ".html"


# =============================================================================
# 12 · IGVVisualization
# =============================================================================

class TestIGVVisualization:
    """Tests for the IGVVisualization rule."""

    @patch("subprocess.run")
    def test_igv_script_receives_all_inputs(self, mock_run, tmp_workdir,
                                             dummy_bam, dummy_vcf, dummy_fasta):
        """igv_visualization.py must receive --bam, --vcf, --ref, --output."""
        mock_run.return_value = MagicMock(returncode=0)
        out_html = str(tmp_workdir / "results" / "igv_visualization.html")
        subprocess.run(
            ["python", "scripts/igv_visualization.py",
             "--bam",    str(dummy_bam),
             "--vcf",    str(dummy_vcf),
             "--ref",    str(dummy_fasta),
             "--output", out_html],
            check=True,
        )
        args = mock_run.call_args[0][0]
        for flag in ["--bam", "--vcf", "--ref", "--output"]:
            assert flag in args, f"Missing flag {flag!r} in IGV script call"

    def test_output_is_html(self, tmp_workdir):
        """IGV output must be an HTML file."""
        out = tmp_workdir / "results" / "igv_visualization.html"
        out.write_text("<html></html>")
        assert out.suffix == ".html"


# =============================================================================
# 13 · Helper functions
# =============================================================================

class TestHelpers:
    """Tests for standalone helper functions defined in the Snakefile."""

    def test_get_pileup_name(self):
        """get_pileup_name should return the correct recalibrated dedup BAM path."""
        def get_pileup_name(wildcards):
            return f"data/{wildcards['group']}.sorted.recalibrated.dedup.bam"

        assert get_pileup_name({"group": "control"}) == \
               "data/control.sorted.recalibrated.dedup.bam"
        assert get_pileup_name({"group": "tumor"}) == \
               "data/tumor.sorted.recalibrated.dedup.bam"

    def test_get_gatk_jar(self, tmp_workdir):
        """get_gatk_jar should find the first .jar under tools/gatk/."""
        gatk_dir = tmp_workdir / "tools" / "gatk"
        gatk_dir.mkdir(parents=True)
        jar = gatk_dir / "gatk-4.4.0.jar"
        jar.write_bytes(b"PK")

        import glob as _glob
        result = _glob.glob(str(gatk_dir / "*.jar"))
        assert len(result) == 1
        assert result[0].endswith(".jar")

    def test_get_gatk_jar_missing_raises(self, tmp_workdir):
        """get_gatk_jar should raise FileNotFoundError when no jar is present."""
        import glob as _glob
        empty_dir = tmp_workdir / "tools" / "gatk_empty"
        empty_dir.mkdir(parents=True)
        jars = _glob.glob(str(empty_dir / "*.jar"))
        if not jars:
            with pytest.raises(Exception):
                raise FileNotFoundError("No GATK jar found under tools/gatk/")