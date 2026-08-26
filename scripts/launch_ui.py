#!/usr/bin/env python3
"""
=============================================================================
launch_ui.py (ONT Human Variation Suite)
Standalone, zero-dependency Python 3 Web UI server for configuring global sample
paths, reference paths, and parameters, validating directories with 'ls -l'
sneak-peek, native OS file/folder browsing, auto-detecting Epi2ME callsets,
and running single or multi-pipeline suites directly from the browser.

Usage:
  python3 scripts/launch_ui.py
  python3 scripts/launch_ui.py --port 8080 --no-browser
=============================================================================
"""

import os
import sys
import json
import glob
import re
import traceback
import argparse
import webbrowser
import subprocess
import threading
from http.server import HTTPServer, ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

CURRENT_PROCESS = None
CURRENT_PROCESS_LOCK = threading.Lock()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
HTML_FILE = os.path.join(SCRIPT_DIR, "ui_index.html")

PIPELINES = ["sv", "methylation", "cnv", "snv"]

def sanitize_prefix(prefix):
    """Clean raw sample prefix by removing trailing pipeline extensions if pasted accidentally."""
    if not prefix:
        return ""
    prefix = prefix.strip()
    # Strip common file suffixes if user pasted a full filename
    for ext in [".wf_sv.vcf.gz", ".wf_mods.bedmethyl.gz", ".wf_cnv.vcf.gz", ".wf_snp.vcf.gz", ".wf_snp_clinvar.vcf.gz", ".wf_", ".vcf.gz", ".vcf", ".gz"]:
        if prefix.endswith(ext):
            prefix = prefix[:-len(ext)]
    return prefix.rstrip("._")

def auto_detect_prefix(input_dir):
    """Scan input_dir for Epi2ME files and extract the true raw sample prefix."""
    if not input_dir or not os.path.exists(input_dir):
        return ""
    
    # Search for *.wf_*.vcf.gz or *.bedmethyl.gz
    patterns = [
        "*.wf_sv.vcf.gz",
        "*.wf_mods.bedmethyl.gz",
        "*.wf_cnv.vcf.gz",
        "*.wf_snp.vcf.gz"
    ]
    for pat in patterns:
        matches = glob.glob(os.path.join(input_dir, pat))
        if matches:
            basename = os.path.basename(matches[0])
            return sanitize_prefix(basename)
    return ""

def clean_base_path(path_str, sample_id):
    """Strip compounding subfolders like /<sample_id>/<pipe> from base output/work dir."""
    if not path_str:
        return ""
    cleaned = path_str.strip().rstrip("/")
    for p in PIPELINES + ["results", "work"]:
        pattern = rf"(/{re.escape(sample_id)})?(/{p})+$"
        cleaned = re.sub(pattern, "", cleaned)
        if cleaned.endswith(f"/{sample_id}"):
            cleaned = cleaned[:-len(sample_id)-1]
    return cleaned.rstrip("/")

def check_tool(name):
    import shutil
    path = shutil.which(name)
    return {"installed": path is not None, "path": path or ""}

class ReusableHTTPServer(ThreadingHTTPServer):
    allow_reuse_address = True
    daemon_threads = True

class SuiteUIHandler(BaseHTTPRequestHandler):

    def _send_json(self, data, status=200):
        try:
            body = json.dumps(data).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        except Exception as e:
            sys.stderr.write(f"Error sending JSON response: {e}\n")

    def _send_html(self, content, status=200):
        try:
            body = content.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        except Exception as e:
            sys.stderr.write(f"Error sending HTML response: {e}\n")

    def do_OPTIONS(self):
        """Handle CORS preflight requests from browser fetch()."""
        self.send_response(200, "ok")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def do_GET(self):
        try:
            parsed = urlparse(self.path)
            path = parsed.path

            if path == "/" or path == "/index.html":
                if os.path.exists(HTML_FILE):
                    with open(HTML_FILE, "r", encoding="utf-8") as f:
                        self._send_html(f.read())
                else:
                    self._send_html("<h1>Error: ui_index.html not found</h1>", 404)

            elif path == "/api/config":
                baseline_cfg = os.path.join(REPO_DIR, "pipelines", "sv", "config", "pipeline_config.yaml")
                baseline_ref = os.path.join(REPO_DIR, "pipelines", "sv", "config", "reference_paths.yaml")

                if not os.path.exists(baseline_cfg):
                    baseline_cfg = os.path.join(REPO_DIR, "pipelines", "sv", "config", "pipeline_config.example.yaml")
                if not os.path.exists(baseline_ref):
                    baseline_ref = os.path.join(REPO_DIR, "pipelines", "sv", "config", "reference_paths.example.yaml")

                cfg_data = {}
                if os.path.exists(baseline_cfg):
                    cfg_data.update(self._parse_simple_yaml(baseline_cfg))
                if os.path.exists(baseline_ref):
                    cfg_data.update(self._parse_simple_yaml(baseline_ref))

                self._send_json({"ok": True, "config": cfg_data})

            elif path == "/api/load-fixtures":
                fixture_dir = os.path.join(REPO_DIR, "pipelines", "sv", "tests", "fixtures")
                gene_bed = os.path.join(fixture_dir, "dummy_genes.bed")
                out_dir = os.path.join(REPO_DIR, "results", "suite_sample_01")
                wrk_dir = os.path.join(REPO_DIR, "work", "suite_sample_01")

                self._send_json({
                    "ok": True,
                    "config": {
                        "sample_id": "TEST_SAMPLE",
                        "raw_sample_prefix": "test_sample",
                        "run_name": "Test Fixtures Verification Run",
                        "input_dir": fixture_dir,
                        "output_dir": out_dir,
                        "work_dir": wrk_dir,
                        "gene_bed": gene_bed,
                        "archive_root": ""
                    }
                })

            else:
                self._send_json({"error": "Not Found"}, 404)

        except Exception as e:
            sys.stderr.write(f"GET error: {traceback.format_exc()}\n")
            self._send_json({"ok": False, "error": str(e)}, 500)

    def do_POST(self):
        global CURRENT_PROCESS
        try:
            parsed = urlparse(self.path)
            length_hdr = self.headers.get("Content-Length")
            length = int(length_hdr) if length_hdr else 0
            body = self.rfile.read(length).decode("utf-8") if length > 0 else "{}"
            try:
                payload = json.loads(body)
            except Exception:
                payload = {}

            if parsed.path == "/api/browse-path":
                browse_type = payload.get("type", "folder")
                prompt = payload.get("prompt", "Select Location")
                selected_path = ""

                if sys.platform == "darwin":
                    try:
                        if browse_type == "file":
                            cmd = f'osascript -e \'POSIX path of (choose file with prompt "{prompt}")\''
                        else:
                            cmd = f'osascript -e \'POSIX path of (choose folder with prompt "{prompt}")\''
                        res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                        if res.returncode == 0:
                            selected_path = res.stdout.strip()
                    except Exception as e:
                        sys.stderr.write(f"osascript error: {e}\n")

                if not selected_path and threading.current_thread() is threading.main_thread():
                    try:
                        import tkinter as tk
                        from tkinter import filedialog
                        root = tk.Tk()
                        root.withdraw()
                        root.attributes('-topmost', True)
                        if browse_type == "file":
                            selected_path = filedialog.askopenfilename(title=prompt)
                        else:
                            selected_path = filedialog.askdirectory(title=prompt)
                        root.destroy()
                    except Exception:
                        pass

                self._send_json({"ok": True, "path": selected_path})

            elif parsed.path == "/api/save-config":
                selected_pipelines = payload.get("pipelines", ["sv"])
                sample_id = payload.get("sample_id", "SAMPLE_01").strip()
                raw_prefix = payload.get("raw_sample_prefix", "").strip()
                run_name = payload.get("run_name", "Epi2ME Downstream Analysis").strip()
                input_dir = payload.get("input_dir", "").strip().rstrip("/")
                raw_out_dir = payload.get("output_dir", "").strip()
                raw_wrk_dir = payload.get("work_dir", "").strip()
                gene_bed = payload.get("gene_bed", "").strip()
                archive_root = payload.get("archive_root", "").strip().rstrip("/")

                # Prevent path compounding
                base_output_dir = clean_base_path(raw_out_dir, sample_id)
                base_work_dir = clean_base_path(raw_wrk_dir, sample_id)

                # Auto-detect true raw sample prefix from input_dir if possible
                detected_prefix = auto_detect_prefix(input_dir) if input_dir else ""
                if detected_prefix:
                    raw_prefix = detected_prefix
                else:
                    raw_prefix = sanitize_prefix(raw_prefix) or sample_id

                # Resolve absolute gene_bed path if provided
                if gene_bed and os.path.exists(gene_bed):
                    gene_bed = os.path.abspath(gene_bed)

                min_qual = payload.get("min_qual", "20")
                min_vaf = payload.get("min_vaf", "0.15")
                min_dp = payload.get("min_dp", "10")
                pass_only = payload.get("pass_only", "true")

                # Resolve FASTA reference for SV if available
                sv_fasta = ""
                if input_dir and os.path.exists(input_dir):
                    fasta_candidates = glob.glob(os.path.join(input_dir, "*.fa")) + glob.glob(os.path.join(input_dir, "*.fasta")) + glob.glob(os.path.join(input_dir, "*.fna"))
                    if fasta_candidates:
                        sv_fasta = fasta_candidates[0]
                if not sv_fasta:
                    fixture_fa = os.path.join(REPO_DIR, "pipelines", "sv", "tests", "fixtures", "dummy_reference.fa")
                    if os.path.exists(fixture_fa):
                        sv_fasta = fixture_fa

                # Resolve Methylation reference files
                meth_fixture_dir = os.path.join(REPO_DIR, "pipelines", "methylation", "tests", "fixtures")
                chrom_sizes = os.path.join(meth_fixture_dir, "test_chrom.sizes")
                promoter_bed = os.path.join(meth_fixture_dir, "test_promoters.bed")
                cpg_island_bed = os.path.join(meth_fixture_dir, "test_cpg_islands.bed")

                saved_files = []
                for pipe in selected_pipelines:
                    if pipe not in PIPELINES:
                        continue
                    pipe_dir = os.path.join(REPO_DIR, "pipelines", pipe)
                    cfg_dir = os.path.join(pipe_dir, "config")
                    os.makedirs(cfg_dir, exist_ok=True)

                    cfg_path = os.path.join(cfg_dir, "pipeline_config.yaml")
                    ref_path = os.path.join(cfg_dir, "reference_paths.yaml")

                    # Output directory clean hierarchy: <base_output_dir>/<sample_id>/<pipe>
                    if base_output_dir:
                        out_dir = os.path.join(base_output_dir, sample_id, pipe)
                        wrk_dir = os.path.join(base_work_dir, sample_id, pipe) if base_work_dir else os.path.join(pipe_dir, "work")
                    else:
                        out_dir = os.path.join(pipe_dir, "results")
                        wrk_dir = os.path.join(pipe_dir, "work")

                    with open(cfg_path, "w", encoding="utf-8") as f:
                        f.write(f"# Pipeline Config for {pipe.upper()}\n")
                        f.write(f"sample:\n")
                        f.write(f'  sample_id: "{sample_id}"\n')
                        f.write(f'  raw_sample_prefix: "{raw_prefix}"\n')
                        f.write(f'  run_description: "{run_name}"\n\n')
                        f.write(f"paths:\n")
                        f.write(f'  input_dir: "{input_dir}"\n')
                        f.write(f'  output_dir: "{out_dir}"\n')
                        f.write(f'  work_dir: "{wrk_dir}"\n')
                        f.write(f'  repo_dir: "{pipe_dir}"\n\n')
                        f.write(f"reference:\n")
                        f.write(f'  config_file: "config/reference_paths.yaml"\n\n')
                        f.write(f"input_files:\n")
                        if pipe == "sv":
                            f.write(f'  sv_vcf: "${{sample.raw_sample_prefix}}.wf_sv.vcf.gz"\n')
                            f.write(f'  cnv_vcf: "${{sample.raw_sample_prefix}}.wf_cnv.vcf.gz"\n')
                            f.write(f'  snv_vcf: "${{sample.raw_sample_prefix}}.wf_snp.vcf.gz"\n')
                            f.write(f'  snv_vcf_clinvar: "${{sample.raw_sample_prefix}}.wf_snp_clinvar.vcf.gz"\n')
                            f.write(f'  str_vcf: "${{sample.raw_sample_prefix}}.wf_str.vcf.gz"\n\n')
                            f.write(f"filtering:\n")
                            f.write(f'  min_sv_length: 30\n')
                            f.write(f'  min_read_support: 4\n')
                            f.write(f'  min_qual: {min_qual}\n')
                            f.write(f'  min_vaf: {min_vaf}\n')
                            f.write(f'  pass_only: {pass_only}\n\n')
                            f.write(f"sv_categories:\n  deletions: [\"DEL\"]\n  insertions: [\"INS\"]\n  duplications: [\"DUP\"]\n  inversions: [\"INV\"]\n  translocations: [\"BND\", \"TRA\"]\n  complex_rearrangements: [\"CPX\", \"BND_CLUSTER\"]\n\n")
                        elif pipe == "methylation":
                            f.write(f'  mod_bed: "${{sample.raw_sample_prefix}}.wf_mods.bedmethyl.gz"\n')
                            f.write(f'  mod_bedmethyl: "${{sample.raw_sample_prefix}}.wf_mods.bedmethyl.gz"\n')
                            f.write(f'  mod_bedmethyl_hap1: "${{sample.raw_sample_prefix}}.wf_mods.1.bedmethyl.gz"\n')
                            f.write(f'  mod_bedmethyl_hap2: "${{sample.raw_sample_prefix}}.wf_mods.2.bedmethyl.gz"\n')
                            f.write(f'  mod_bedmethyl_ungrouped: "${{sample.raw_sample_prefix}}.wf_mods.ungrouped.bedmethyl.gz"\n\n')
                            f.write(f"modifications:\n  phased: false\n  primary_mod_code: \"m\"\n  expected_columns: 18\n\n")
                            f.write(f"filtering:\n  min_coverage: {min_dp}\n  max_coverage: 0\n  primary_contigs_only: true\n  include_contigs_regex: \"^chr([1-9]|1[0-9]|2[0-2]|X|Y)$\"\n  expected_primary_contigs: 24\n\n")
                            f.write(f"methylation_states:\n  unmethylated_max_percent: 20\n  methylated_min_percent: 80\n\n")
                            f.write(f"annotation:\n  min_cpgs_per_feature: 5\n  gene_key: \"gene_id\"\n  keep_annotated_cpgs: false\n\n")
                        elif pipe == "cnv":
                            f.write(f'  cnv_vcf: "${{sample.raw_sample_prefix}}.wf_cnv.vcf.gz"\n\n')
                            f.write(f"filtering:\n  min_cnv_length: 10000\n  min_qual: {min_qual}\n  pass_only: {pass_only}\n  max_cn_loss: 1\n  min_cn_gain: 3\n\n")
                            f.write(f"cnv_categories:\n  deletions: [\"DEL\", \"LOSS\"]\n  duplications: [\"DUP\", \"GAIN\", \"AMP\"]\n\n")
                        elif pipe == "snv":
                            f.write(f'  snv_vcf: "${{sample.raw_sample_prefix}}.wf_snp.vcf.gz"\n')
                            f.write(f'  snv_vcf_clinvar: "${{sample.raw_sample_prefix}}.wf_snp_clinvar.vcf.gz"\n\n')
                            f.write(f"filtering:\n  min_dp: {min_dp}\n  min_qual: {min_qual}\n  min_vaf: {min_vaf}\n  pass_only: {pass_only}\n\n")
                            f.write(f"snv_categories:\n  snvs: [\"SNV\", \"SNP\"]\n  indels: [\"INDEL\", \"INS\", \"DEL\"]\n\n")
                            f.write(f"clinvar:\n  extract_pathogenic: true\n\n")

                        f.write(f"archive:\n")
                        f.write(f'  archive_root: "{archive_root}"\n')
                        f.write(f'  compress: false\n\n')
                        f.write(f"compute:\n  threads: 8\n\nlogging:\n  log_dir: \"logs\"\n")

                    with open(ref_path, "w", encoding="utf-8") as f:
                        if pipe == "sv":
                            f.write(f"genome:\n  fasta: \"{sv_fasta}\"\n  build: \"GRCh38\"\n\nannotation:\n")
                            f.write(f'  gene_bed: "{gene_bed}"\n\n')
                            f.write(f"tools:\n  bcftools: \"bcftools\"\n  bedtools: \"bedtools\"\n  tabix: \"tabix\"\n  R: \"Rscript\"\n")
                        elif pipe == "methylation":
                            f.write(f"genome:\n  fasta: \"\"\n  build: \"GRCh38\"\n  chrom_sizes: \"{chrom_sizes}\"\n\nannotation:\n")
                            f.write(f'  gene_bed: "{gene_bed}"\n')
                            f.write(f'  promoter_bed: "{promoter_bed}"\n')
                            f.write(f'  cpg_island_bed: "{cpg_island_bed}"\n\n')
                            f.write(f"tools:\n  bedtools: \"bedtools\"\n  tabix: \"tabix\"\n  bgzip: \"bgzip\"\n  R: \"Rscript\"\n")
                        else:
                            f.write(f"genome:\n  build: \"GRCh38\"\n\nannotation:\n")
                            f.write(f'  gene_bed: "{gene_bed}"\n\n')
                            f.write(f"tools:\n  bcftools: \"bcftools\"\n  bedtools: \"bedtools\"\n  tabix: \"tabix\"\n  R: \"Rscript\"\n")

                    saved_files.append(pipe)

                self._send_json({
                    "ok": True,
                    "saved_pipelines": saved_files,
                    "sanitized_raw_prefix": raw_prefix,
                    "clean_output_dir": base_output_dir,
                    "clean_work_dir": base_work_dir
                })

            elif parsed.path == "/api/validate-paths":
                input_dir = payload.get("input_dir", "").strip().rstrip("/")
                gene_bed = payload.get("gene_bed", "").strip()
                output_dir = payload.get("output_dir", "").strip().rstrip("/")
                work_dir = payload.get("work_dir", "").strip().rstrip("/")
                raw_prefix = payload.get("raw_sample_prefix", "").strip()
                selected_pipelines = payload.get("pipelines", PIPELINES)

                detected_prefix = auto_detect_prefix(input_dir) if input_dir else ""
                clean_prefix = detected_prefix or sanitize_prefix(raw_prefix)

                out_parent = os.path.dirname(output_dir) or "."
                wrk_parent = os.path.dirname(work_dir) or "."

                validation = {
                    "input_dir": {
                        "path": input_dir,
                        "exists": os.path.isdir(input_dir) if input_dir else False,
                        "readable": os.access(input_dir, os.R_OK) if input_dir and os.path.exists(input_dir) else False,
                        "detected_prefix": detected_prefix
                    },
                    "callsets": {},
                    "gene_bed": {
                        "path": gene_bed,
                        "exists": os.path.isfile(gene_bed) if gene_bed else False,
                        "readable": os.access(gene_bed, os.R_OK) if gene_bed and os.path.exists(gene_bed) else False
                    },
                    "output_dir": {
                        "path": output_dir,
                        "exists": os.path.exists(output_dir) if output_dir else False,
                        "writable": os.access(output_dir, os.W_OK) if output_dir and os.path.exists(output_dir) else (os.access(out_parent, os.W_OK) if output_dir else False)
                    },
                    "work_dir": {
                        "path": work_dir,
                        "exists": os.path.exists(work_dir) if work_dir else False,
                        "writable": os.access(work_dir, os.W_OK) if work_dir and os.path.exists(work_dir) else (os.access(wrk_parent, os.W_OK) if work_dir else False)
                    },
                    "tools": {
                        "bcftools": check_tool("bcftools"),
                        "bedtools": check_tool("bedtools"),
                        "tabix": check_tool("tabix"),
                        "Rscript": check_tool("Rscript")
                    }
                }

                expected_filenames = {
                    "sv": f"{clean_prefix}.wf_sv.vcf.gz",
                    "methylation": f"{clean_prefix}.wf_mods.bedmethyl.gz",
                    "cnv": f"{clean_prefix}.wf_cnv.vcf.gz",
                    "snv": f"{clean_prefix}.wf_snp.vcf.gz"
                }

                for pipe in PIPELINES:
                    fname = expected_filenames.get(pipe, "")
                    fpath = os.path.join(input_dir, fname) if (input_dir and fname) else ""
                    exists = os.path.isfile(fpath) if fpath else False
                    size = os.path.getsize(fpath) if exists else 0
                    validation["callsets"][pipe] = {
                        "expected_file": fname,
                        "path": fpath,
                        "exists": exists,
                        "size_bytes": size,
                        "selected": pipe in selected_pipelines
                    }

                self._send_json({"ok": True, "validation": validation})

            elif parsed.path == "/api/peek-path":
                input_dir = payload.get("input_dir", "").strip().rstrip("/")
                gene_bed = payload.get("gene_bed", "").strip()
                output_dir = payload.get("output_dir", "").strip().rstrip("/")
                raw_prefix = payload.get("raw_sample_prefix", "").strip()

                detected_prefix = auto_detect_prefix(input_dir) if input_dir else ""
                clean_prefix = detected_prefix or sanitize_prefix(raw_prefix)

                results = {}
                if input_dir and os.path.exists(input_dir):
                    try:
                        ls_res = subprocess.run(["ls", "-lh", input_dir], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                        dir_contents = ls_res.stdout
                        matching = glob.glob(os.path.join(input_dir, "*.wf_*.gz")) + glob.glob(os.path.join(input_dir, "*.bedmethyl.gz"))
                        matching_names = [os.path.basename(m) for m in matching]

                        results["input_dir"] = {
                            "exists": True,
                            "ls_output": dir_contents[:2000],
                            "matching_files": matching_names,
                            "detected_prefix": detected_prefix,
                            "sanitized_prefix": clean_prefix
                        }
                    except Exception as e:
                        results["input_dir"] = {"exists": True, "error": str(e)}
                else:
                    results["input_dir"] = {"exists": False}

                if gene_bed:
                    results["gene_bed"] = {"exists": os.path.exists(gene_bed)}
                if output_dir:
                    results["output_dir"] = {"exists": os.path.exists(output_dir)}

                self._send_json({"ok": True, "peek": results})

            elif parsed.path == "/api/check-deps":
                selected_pipelines = payload.get("pipelines", ["sv"])
                outputs = []
                for pipe in selected_pipelines:
                    cmd = ["bash", os.path.join(REPO_DIR, "scripts", "check_dependencies.sh"), "--pipeline", pipe]
                    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                    outputs.append(f"=== Dependency Check [{pipe.upper()}] ===\n{res.stdout}")
                self._send_json({"ok": True, "output": "\n\n".join(outputs)})

            elif parsed.path == "/api/stop-suite":
                stopped = False
                with CURRENT_PROCESS_LOCK:
                    if CURRENT_PROCESS and CURRENT_PROCESS.poll() is None:
                        try:
                            CURRENT_PROCESS.terminate()
                            try:
                                CURRENT_PROCESS.wait(timeout=3)
                            except subprocess.TimeoutExpired:
                                CURRENT_PROCESS.kill()
                            stopped = True
                        except Exception as e:
                            sys.stderr.write(f"Error stopping process: {e}\n")
                        CURRENT_PROCESS = None
                self._send_json({"ok": True, "stopped": stopped})

            elif parsed.path == "/api/run-suite-stream":
                selected_pipelines = payload.get("pipelines", ["sv"])
                
                self.send_response(200)
                self.send_header("Content-Type", "text/event-stream; charset=utf-8")
                self.send_header("Cache-Control", "no-cache")
                self.send_header("Access-Control-Allow-Origin", "*")
                self.send_header("Connection", "keep-alive")
                self.end_headers()

                def send_event(event_type, msg, **kwargs):
                    try:
                        data_dict = {"type": event_type, "message": msg}
                        data_dict.update(kwargs)
                        data = json.dumps(data_dict)
                        self.wfile.write(f"data: {data}\n\n".encode("utf-8"))
                        self.wfile.flush()
                    except Exception:
                        pass

                send_event("start", f"Starting suite execution for: {', '.join(selected_pipelines)}", pipelines=selected_pipelines)

                total_pipes = len(selected_pipelines)
                for p_idx, pipe in enumerate(selected_pipelines):
                    pipe_dir = os.path.join(REPO_DIR, "pipelines", pipe)
                    script = os.path.join(pipe_dir, "scripts", "bash", "04_run_all.sh")
                    cfg = os.path.join(pipe_dir, "config", "pipeline_config.yaml")

                    if not os.path.exists(cfg):
                        send_event("error", f"❌ [{pipe.upper()}] Config not found: {cfg}. Save config first.", pipeline=pipe)
                        continue

                    send_event("pipeline_start", f"\n====================================================\n▶ LAUNCHING PIPELINE: {pipe.upper()}\n====================================================", pipeline=pipe, pipe_index=p_idx, total_pipelines=total_pipes)
                    
                    cmd = ["bash", script, cfg]
                    proc = subprocess.Popen(cmd, cwd=pipe_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

                    with CURRENT_PROCESS_LOCK:
                        CURRENT_PROCESS = proc

                    for line in iter(proc.stdout.readline, ""):
                        if not line:
                            break
                        line_str = line.rstrip()

                        stage_match = re.search(r'>>> Stage (\w+):\s*(.*)', line_str)
                        if stage_match:
                            stage_id = stage_match.group(1)
                            stage_desc = stage_match.group(2)
                            send_event("stage_start", line_str, pipeline=pipe, stage_id=stage_id, stage_desc=stage_desc, pipe_index=p_idx, total_pipelines=total_pipes)
                        else:
                            send_event("log", line_str, pipeline=pipe)

                    proc.stdout.close()
                    rc = proc.wait()

                    with CURRENT_PROCESS_LOCK:
                        CURRENT_PROCESS = None

                    if rc != 0:
                        send_event("pipeline_error", f"❌ [{pipe.upper()}] Pipeline failed with exit code {rc}. Aborting suite run.", pipeline=pipe, rc=rc)
                        send_event("end", f"Suite execution failed at pipeline: {pipe}", success=False)
                        return

                    send_event("pipeline_complete", f"✓ [{pipe.upper()}] Pipeline completed successfully.", pipeline=pipe, pipe_index=p_idx, total_pipelines=total_pipes)

                send_event("end", "✓ Full pipeline suite completed successfully!", success=True)
                return

            elif parsed.path == "/api/run-suite":
                selected_pipelines = payload.get("pipelines", ["sv"])
                suite_logs = []

                for pipe in selected_pipelines:
                    pipe_dir = os.path.join(REPO_DIR, "pipelines", pipe)
                    script = os.path.join(pipe_dir, "scripts", "bash", "04_run_all.sh")
                    cfg = os.path.join(pipe_dir, "config", "pipeline_config.yaml")

                    if not os.path.exists(cfg):
                        suite_logs.append(f"❌ [{pipe.upper()}] Config not found: {cfg}. Save config first.")
                        continue

                    suite_logs.append(f"====================================================\n▶ LAUNCHING PIPELINE: {pipe.upper()}\n====================================================")
                    cmd = ["bash", script, cfg]
                    res = subprocess.run(cmd, cwd=pipe_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                    suite_logs.append(res.stdout)
                    if res.returncode != 0:
                        suite_logs.append(f"❌ [{pipe.upper()}] Pipeline failed with exit code {res.returncode}. Aborting suite run.")
                        break

                self._send_json({"ok": True, "output": "\n".join(suite_logs)})

            else:
                self._send_json({"ok": False, "error": "Not Found"}, 404)

        except Exception as e:
            sys.stderr.write(f"POST error: {traceback.format_exc()}\n")
            self._send_json({"ok": False, "error": str(e)}, 500)

    def _parse_simple_yaml(self, filepath):
        result = {}
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#") or ":" not in line:
                        continue
                    parts = line.split(":", 1)
                    k = parts[0].strip()
                    v = parts[1].strip().strip('"').strip("'")
                    if k and v:
                        result[k] = v
        except Exception:
            pass
        return result

    def log_message(self, format, *args):
        return

def main():
    parser = argparse.ArgumentParser(description="Launch ONT Human Variation Suite Web UI")
    parser.add_argument("--port", type=int, default=5005)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()

    httpd = None
    candidate_ports = [args.port, 5005, 5006, 5007, 8080, 8081, 8888]
    seen = set()
    unique_ports = [p for p in candidate_ports if not (p in seen or seen.add(p))]

    for p in unique_ports:
        for host in ["127.0.0.1", "localhost", ""]:
            try:
                server_address = (host, p)
                httpd = ReusableHTTPServer(server_address, SuiteUIHandler)
                port = p
                break
            except OSError:
                continue
        if httpd:
            break

    if not httpd:
        print("Error: Could not bind to any available port.")
        sys.exit(1)

    url = f"http://localhost:{port}"
    print("============================================================")
    print("  ONT Human Variation Suite — Global Web UI Control")
    print("============================================================")
    print(f"  Server listening on: {url}")
    print("  Press Ctrl+C to stop the server.")
    print("============================================================")

    if not args.no_browser:
        try:
            webbrowser.open(url)
        except Exception:
            pass

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping Web UI server...")
        httpd.server_close()

if __name__ == "__main__":
    main()
