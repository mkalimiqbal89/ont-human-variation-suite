#!/usr/bin/env python3
"""
=============================================================================
launch_ui.py (ONT Human Variation Suite)
Standalone, zero-dependency Python 3 Web UI server for configuring sample paths,
reference paths, and parameters, and running pipelines directly from the browser.

Usage:
  python3 scripts/launch_ui.py
  python3 scripts/launch_ui.py --pipeline cnv
  python3 scripts/launch_ui.py --port 8080 --no-browser
=============================================================================
"""

import os
import sys
import json
import argparse
import webbrowser
import subprocess
from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
HTML_FILE = os.path.join(SCRIPT_DIR, "ui_index.html")

class SuiteUIHandler(BaseHTTPRequestHandler):

    def _send_json(self, data, status=200):
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(json.dumps(data).encode("utf-8"))

    def _send_html(self, content, status=200):
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.end_headers()
        self.wfile.write(content.encode("utf-8"))

    def do_GET(self):
        parsed = urlparse(self.path)
        path = parsed.path
        query = parse_qs(parsed.query)

        if path == "/" or path == "/index.html":
            if os.path.exists(HTML_FILE):
                with open(HTML_FILE, "r", encoding="utf-8") as f:
                    self._send_html(f.read())
            else:
                self._send_html("<h1>Error: ui_index.html not found</h1>", 404)

        elif path == "/api/config":
            pipeline = query.get("pipeline", ["sv"])[0]
            pipe_dir = os.path.join(REPO_DIR, "pipelines", pipeline)
            cfg_file = os.path.join(pipe_dir, "config", "pipeline_config.yaml")
            ex_cfg_file = os.path.join(pipe_dir, "config", "pipeline_config.example.yaml")
            ref_cfg_file = os.path.join(pipe_dir, "config", "reference_paths.yaml")
            ex_ref_cfg_file = os.path.join(pipe_dir, "config", "reference_paths.example.yaml")

            target_cfg = cfg_file if os.path.exists(cfg_file) else ex_cfg_file
            target_ref = ref_cfg_file if os.path.exists(ref_cfg_file) else ex_ref_cfg_file

            cfg_data = {}
            if os.path.exists(target_cfg):
                cfg_data.update(self._parse_yaml(target_cfg))
            if os.path.exists(target_ref):
                cfg_data.update(self._parse_yaml(target_ref))

            self._send_json({"ok": True, "config": cfg_data})

        else:
            self._send_json({"error": "Not Found"}, 404)

    def do_POST(self):
        parsed = urlparse(self.path)
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length).decode("utf-8") if length > 0 else "{}"
        try:
            payload = json.loads(body)
        except Exception:
            payload = {}

        if parsed.path == "/api/config":
            pipeline = payload.get("pipeline", "sv")
            pipe_dir = os.path.join(REPO_DIR, "pipelines", pipeline)
            cfg_dir = os.path.join(pipe_dir, "config")
            os.makedirs(cfg_dir, exist_ok=True)

            cfg_path = os.path.join(cfg_dir, "pipeline_config.yaml")
            ref_path = os.path.join(cfg_dir, "reference_paths.yaml")

            # Write pipeline_config.yaml
            with open(cfg_path, "w", encoding="utf-8") as f:
                f.write(f"sample:\n")
                f.write(f'  sample_id: "{payload.get("sample_id", "SAMPLE_01")}"\n')
                f.write(f'  raw_sample_prefix: "{payload.get("raw_sample_prefix", "sample_01")}"\n\n')
                f.write(f"paths:\n")
                f.write(f'  input_dir: "{payload.get("input_dir", "")}"\n')
                f.write(f'  output_dir: "{payload.get("output_dir", "")}"\n')
                f.write(f'  work_dir: "{payload.get("work_dir", "")}"\n')
                f.write(f'  repo_dir: "{pipe_dir}"\n\n')
                f.write(f"reference:\n")
                f.write(f'  config_file: "config/reference_paths.yaml"\n\n')
                f.write(f"input_files:\n")
                if pipeline == "sv":
                    f.write(f'  sv_vcf: "${{sample.raw_sample_prefix}}.wf_sv.vcf.gz"\n')
                elif pipeline == "methylation":
                    f.write(f'  mod_bed: "${{sample.raw_sample_prefix}}.wf_mods.bedmethyl.gz"\n')
                elif pipeline == "cnv":
                    f.write(f'  cnv_vcf: "${{sample.raw_sample_prefix}}.wf_cnv.vcf.gz"\n')
                elif pipeline == "snv":
                    f.write(f'  snv_vcf: "${{sample.raw_sample_prefix}}.wf_snp.vcf.gz"\n')
                    f.write(f'  snv_vcf_clinvar: "${{sample.raw_sample_prefix}}.wf_snp_clinvar.vcf.gz"\n')
                f.write(f"\nfiltering:\n")
                f.write(f'  min_qual: {payload.get("min_qual", 20)}\n')
                f.write(f'  min_vaf: {payload.get("min_vaf", 0.15)}\n')
                f.write(f'  min_dp: {payload.get("min_dp", 10)}\n')
                f.write(f'  pass_only: {payload.get("pass_only", "true")}\n\n')
                f.write(f"compute:\n  threads: 8\n\nlogging:\n  log_dir: \"logs\"\n")

            # Write reference_paths.yaml
            with open(ref_path, "w", encoding="utf-8") as f:
                f.write(f"genome:\n  build: \"GRCh38\"\n\nannotation:\n")
                f.write(f'  gene_bed: "{payload.get("gene_bed", "")}"\n\n')
                f.write(f"tools:\n  bcftools: \"bcftools\"\n  bedtools: \"bedtools\"\n  tabix: \"tabix\"\n  R: \"Rscript\"\n")

            self._send_json({"ok": True, "message": f"Saved config to {cfg_path}"})

        elif parsed.path == "/api/validate-paths":
            check_path = payload.get("path", "")
            exists = os.path.exists(check_path) if check_path else False
            self._send_json({"ok": True, "path": check_path, "exists": exists})

        elif parsed.path == "/api/check-deps":
            pipeline = payload.get("pipeline", "sv")
            cmd = ["bash", os.path.join(REPO_DIR, "scripts", "check_dependencies.sh"), "--pipeline", pipeline]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            self._send_json({"ok": res.returncode == 0, "output": res.stdout})

        elif parsed.path == "/api/run-pipeline":
            pipeline = payload.get("pipeline", "sv")
            pipe_dir = os.path.join(REPO_DIR, "pipelines", pipeline)
            script = os.path.join(pipe_dir, "scripts", "bash", "04_run_all.sh")
            cfg = os.path.join(pipe_dir, "config", "pipeline_config.yaml")

            if not os.path.exists(cfg):
                self._send_json({"ok": False, "output": f"Config not found: {cfg}. Save config first."})
                return

            cmd = ["bash", script, cfg]
            res = subprocess.run(cmd, cwd=pipe_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            self._send_json({"ok": res.returncode == 0, "output": res.stdout})

        else:
            self._send_json({"error": "Not Found"}, 404)

    def _parse_yaml(self, filepath):
        """Simple regex-free scalar YAML key-value parser for simple key-value templates."""
        result = {}
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or ":" not in line:
                    continue
                k, v = line.split(":", 1)
                k = k.strip()
                v = v.strip().strip('"').strip("'")
                if k and v and not k.startswith("-"):
                    result[k] = v
        return result

    def log_message(self, format, *args):
        # Quiet standard logging
        return

def main():
    parser = argparse.ArgumentParser(description="Launch ONT Human Variation Suite Web UI")
    parser.add_argument("--pipeline", default="sv", choices=["sv", "methylation", "cnv", "snv"])
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()

    server_address = ("", args.port)
    httpd = HTTPServer(server_address, SuiteUIHandler)

    url = f"http://localhost:{args.port}"
    print("============================================================")
    print("  ONT Human Variation Suite — Interactive Web UI")
    print("============================================================")
    print(f"  Server listening on: {url}")
    print(f"  Selected pipeline  : {args.pipeline}")
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
