#!/usr/bin/env python3
"""
=============================================================================
launch_ui.py (ONT Human Variation Suite)
Standalone, zero-dependency Python 3 Web UI server for configuring global sample
paths, reference paths, and parameters, validating directories with 'ls -l'
sneak-peek, native OS file/folder browsing, and running single or multi-pipeline
suites directly from the browser.

Usage:
  python3 scripts/launch_ui.py
  python3 scripts/launch_ui.py --port 8080 --no-browser
=============================================================================
"""

import os
import sys
import json
import glob
import traceback
import argparse
import webbrowser
import subprocess
from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
HTML_FILE = os.path.join(SCRIPT_DIR, "ui_index.html")

PIPELINES = ["sv", "methylation", "cnv", "snv"]

class ReusableHTTPServer(HTTPServer):
    allow_reuse_address = True

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

            else:
                self._send_json({"error": "Not Found"}, 404)

        except Exception as e:
            sys.stderr.write(f"GET error: {traceback.format_exc()}\n")
            self._send_json({"ok": False, "error": str(e)}, 500)

    def do_POST(self):
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

                if not selected_path:
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
                sample_id = payload.get("sample_id", "SAMPLE_01")
                raw_prefix = payload.get("raw_sample_prefix", "sample_01")
                run_name = payload.get("run_name", "Epi2ME Downstream Analysis")
                input_dir = payload.get("input_dir", "")
                base_output_dir = payload.get("output_dir", "")
                base_work_dir = payload.get("work_dir", "")
                gene_bed = payload.get("gene_bed", "")
                archive_root = payload.get("archive_root", "")

                min_qual = payload.get("min_qual", "20")
                min_vaf = payload.get("min_vaf", "0.15")
                min_dp = payload.get("min_dp", "10")
                pass_only = payload.get("pass_only", "true")

                saved_files = []
                for pipe in selected_pipelines:
                    if pipe not in PIPELINES:
                        continue
                    pipe_dir = os.path.join(REPO_DIR, "pipelines", pipe)
                    cfg_dir = os.path.join(pipe_dir, "config")
                    os.makedirs(cfg_dir, exist_ok=True)

                    cfg_path = os.path.join(cfg_dir, "pipeline_config.yaml")
                    ref_path = os.path.join(cfg_dir, "reference_paths.yaml")

                    if len(selected_pipelines) > 1 and base_output_dir:
                        out_dir = os.path.join(base_output_dir, sample_id, pipe)
                        wrk_dir = os.path.join(base_work_dir, sample_id, pipe) if base_work_dir else os.path.join(pipe_dir, "work")
                    else:
                        out_dir = base_output_dir if base_output_dir else os.path.join(pipe_dir, "results")
                        wrk_dir = base_work_dir if base_work_dir else os.path.join(pipe_dir, "work")

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
                        elif pipe == "methylation":
                            f.write(f'  mod_bed: "${{sample.raw_sample_prefix}}.wf_mods.bedmethyl.gz"\n')
                        elif pipe == "cnv":
                            f.write(f'  cnv_vcf: "${{sample.raw_sample_prefix}}.wf_cnv.vcf.gz"\n')
                        elif pipe == "snv":
                            f.write(f'  snv_vcf: "${{sample.raw_sample_prefix}}.wf_snp.vcf.gz"\n')
                            f.write(f'  snv_vcf_clinvar: "${{sample.raw_sample_prefix}}.wf_snp_clinvar.vcf.gz"\n')

                        f.write(f"\nfiltering:\n")
                        f.write(f'  min_qual: {min_qual}\n')
                        f.write(f'  min_vaf: {min_vaf}\n')
                        f.write(f'  min_dp: {min_dp}\n')
                        f.write(f'  pass_only: {pass_only}\n\n')
                        f.write(f"archive:\n")
                        f.write(f'  archive_root: "{archive_root}"\n')
                        f.write(f'  compress: false\n\n')
                        f.write(f"compute:\n  threads: 8\n\nlogging:\n  log_dir: \"logs\"\n")

                    with open(ref_path, "w", encoding="utf-8") as f:
                        f.write(f"genome:\n  build: \"GRCh38\"\n\nannotation:\n")
                        f.write(f'  gene_bed: "{gene_bed}"\n\n')
                        f.write(f"tools:\n  bcftools: \"bcftools\"\n  bedtools: \"bedtools\"\n  tabix: \"tabix\"\n  R: \"Rscript\"\n")

                    saved_files.append(pipe)

                self._send_json({"ok": True, "saved_pipelines": saved_files})

            elif parsed.path == "/api/peek-path":
                input_dir = payload.get("input_dir", "")
                gene_bed = payload.get("gene_bed", "")
                output_dir = payload.get("output_dir", "")

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
                            "matching_files": matching_names
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
