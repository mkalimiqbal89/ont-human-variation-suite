#!/usr/bin/env bash
# =============================================================================
# check_dependencies.sh
# Reports which external tools and R packages each pipeline needs, whether they
# are present and actually working, and how to install whatever is missing.
#
# Usage:
#   bash scripts/check_dependencies.sh                  check and report only
#   bash scripts/check_dependencies.sh --install        install what is missing
#   bash scripts/check_dependencies.sh --pipeline sv    check one pipeline only
#   bash scripts/check_dependencies.sh --conda          print a conda recipe
#
# Exit status: 0 if everything required is present and working, 1 otherwise.
# Optional tools never affect the exit status.
#
# WHY THIS EXISTS
# The pipelines already check their own tools at stage 00, but only for the
# pipeline being run, and only once you are already trying to run it. Someone
# cloning this repository should be able to find out in one command what they
# need — including the R packages, which are easy to miss because only the SV
# pipeline uses them and only in stages 05-07.
#
# PRESENCE IS NOT ENOUGH. A binary can sit on PATH and still be unusable
# through a broken shared library, so each tool is executed and both its exit
# status and its output are inspected. That distinction is why the SV pipeline
# once reported [OK] for a bcftools that could not run.
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# shellcheck source=../common/lib_common.sh
LIB_COMMON="${REPO_DIR}/common/lib_common.sh"
if [[ ! -f "${LIB_COMMON}" ]]; then
    echo "[FATAL] Cannot find ${LIB_COMMON}" >&2
    exit 1
fi
source "${LIB_COMMON}"

DO_INSTALL=0
ONLY_PIPELINE=""
SHOW_CONDA=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        '#'*)        shift ;;   # zsh passes inline comments through as arguments
        --install)   DO_INSTALL=1; shift ;;
        --pipeline)  ONLY_PIPELINE="$2"; shift 2 ;;
        --conda)     SHOW_CONDA=1; shift ;;
        -h|--help)   sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)           echo "[FATAL] Unknown option: $1" >&2; exit 2 ;;
    esac
done

# --- What each pipeline needs -----------------------------------------------
# Kept here rather than parsed out of the stage scripts: this file is the
# human-readable answer to "what do I need to install", and a parser would break
# the moment a stage script changed shape.
COMMON_TOOLS="bash awk sort gzip"
METHYL_TOOLS="bedtools tabix bgzip Rscript"
SV_TOOLS="bcftools bedtools tabix Rscript"
OPTIONAL_TOOLS="pigz"
SV_R_PACKAGES="yaml ggplot2 scales"
# The methylation pipeline deliberately requires NO R packages — base graphics
# and stats only — so that it cannot fail at the figure stage on a machine where
# the R library path differs from the login shell's.
METHYL_R_PACKAGES=""

# --- Platform and package manager -------------------------------------------
OS="$(uname -s)"
case "${OS}" in
    Darwin) PLATFORM="macOS" ;;
    Linux)  PLATFORM="Linux" ;;
    *)      PLATFORM="${OS}" ;;
esac

PKG_MGR=""
if command -v brew    >/dev/null 2>&1; then PKG_MGR="brew"
elif command -v apt-get >/dev/null 2>&1; then PKG_MGR="apt"
elif command -v dnf     >/dev/null 2>&1; then PKG_MGR="dnf"
elif command -v yum     >/dev/null 2>&1; then PKG_MGR="yum"
fi
HAVE_CONDA=0
command -v conda >/dev/null 2>&1 && HAVE_CONDA=1

# Package names differ per manager. htslib supplies both tabix and bgzip.
pkg_name() {
    local tool="$1" mgr="$2"
    case "${mgr}:${tool}" in
        brew:tabix|brew:bgzip)   echo "htslib" ;;
        brew:Rscript)            echo "r" ;;
        brew:*)                  echo "${tool}" ;;
        apt:bgzip|apt:tabix)     echo "tabix" ;;
        apt:Rscript)             echo "r-base-core" ;;
        apt:*)                   echo "${tool}" ;;
        dnf:bedtools|yum:bedtools) echo "BEDTools" ;;
        dnf:tabix|dnf:bgzip|yum:tabix|yum:bgzip) echo "htslib" ;;
        dnf:Rscript|yum:Rscript) echo "R-core" ;;
        dnf:*|yum:*)             echo "${tool}" ;;
        conda:tabix|conda:bgzip) echo "htslib" ;;
        conda:Rscript)           echo "r-base" ;;
        conda:*)                 echo "${tool}" ;;
        *)                       echo "${tool}" ;;
    esac
}

# --- Tool checking -----------------------------------------------------------
MISSING_REQUIRED=""
MISSING_OPTIONAL=""
BROKEN=""

check_tool() {
    local tool="$1" required="${2:-required}"
    local result status ver

    # The actual "run it and check for a broken shared library" logic lives in
    # common_probe_tool (common/lib_common.sh) — this function only adds the
    # required-vs-optional bookkeeping check_dependencies.sh itself needs.
    result="$(common_probe_tool "${tool}")"
    status="${result%% *}"
    ver="${result#* }"

    case "${status}" in
        MISS)
            printf '  [MISS] %-10s not found on PATH\n' "${tool}"
            if [[ "${required}" == "required" ]]; then
                case " ${MISSING_REQUIRED} " in *" ${tool} "*) ;; *) MISSING_REQUIRED="${MISSING_REQUIRED} ${tool}" ;; esac
            else
                case " ${MISSING_OPTIONAL} " in *" ${tool} "*) ;; *) MISSING_OPTIONAL="${MISSING_OPTIONAL} ${tool}" ;; esac
            fi
            return 1
            ;;
        FAIL)
            printf '  [FAIL] %-10s on PATH but will not run: %s\n' "${tool}" "${ver}"
            case " ${BROKEN} " in *" ${tool} "*) ;; *) BROKEN="${BROKEN} ${tool}" ;; esac
            return 1
            ;;
        *)
            printf '  [OK]   %-10s %s\n' "${tool}" "${ver}"
            return 0
            ;;
    esac
}

MISSING_R_PACKAGES=""

check_r_packages() {
    local pkgs="$1"
    [[ -z "${pkgs}" ]] && { echo "  (none required)"; return 0; }
    if ! command -v Rscript >/dev/null 2>&1; then
        echo "  [SKIP] Rscript not available; cannot check R packages"
        return 1
    fi
    local p
    for p in ${pkgs}; do
        if Rscript -e "quit(status = as.integer(!requireNamespace('${p}', quietly = TRUE)))" >/dev/null 2>&1; then
            printf '  [OK]   %-10s installed\n' "${p}"
        else
            printf '  [MISS] %-10s not installed\n' "${p}"
            case " ${MISSING_R_PACKAGES} " in *" ${p} "*) ;; *) MISSING_R_PACKAGES="${MISSING_R_PACKAGES} ${p}" ;; esac
        fi
    done
}

# =============================================================================
echo "============================================================"
echo " ONT Human Variation Suite — dependency check"
echo "============================================================"
echo "Platform       : ${PLATFORM} ($(uname -m))"
echo "Package manager: ${PKG_MGR:-none detected}"
echo "conda present  : $([[ ${HAVE_CONDA} -eq 1 ]] && echo yes || echo no)"
echo ""

echo "--- Common (all pipelines) ---"
for t in ${COMMON_TOOLS}; do check_tool "${t}"; done

if [[ -z "${ONLY_PIPELINE}" || "${ONLY_PIPELINE}" == "methylation" ]]; then
    echo ""
    echo "--- Methylation pipeline ---"
    for t in ${METHYL_TOOLS}; do check_tool "${t}"; done
    echo "  R packages:"
    check_r_packages "${METHYL_R_PACKAGES}"
fi

if [[ -z "${ONLY_PIPELINE}" || "${ONLY_PIPELINE}" == "sv" ]]; then
    echo ""
    echo "--- SV pipeline ---"
    for t in ${SV_TOOLS}; do check_tool "${t}"; done
    echo "  R packages (used by stages 05-07):"
    check_r_packages "${SV_R_PACKAGES}"
fi

echo ""
echo "--- Optional (performance only) ---"
for t in ${OPTIONAL_TOOLS}; do check_tool "${t}" optional; done
echo "  pigz is a multi-threaded gzip. Without it the pipelines still work;"
echo "  reading a 600 MB+ bedMethyl is simply slower."

# =============================================================================
# Summary and remediation
# =============================================================================
echo ""
echo "============================================================"

NEED_TOOLS="$(echo "${MISSING_REQUIRED} ${BROKEN}" | tr ' ' '\n' | grep -v '^$' | sort -u | tr '\n' ' ')"
NEED_R="$(echo "${MISSING_R_PACKAGES}" | tr ' ' '\n' | grep -v '^$' | sort -u | tr '\n' ' ')"

if [[ -z "${NEED_TOOLS}" && -z "${NEED_R}" ]]; then
    echo " Everything required is present and working."
    echo "============================================================"
    [[ -n "${MISSING_OPTIONAL}" ]] && echo "Optional not installed:${MISSING_OPTIONAL}"
    exit 0
fi

echo " Missing or broken"
echo "============================================================"
[[ -n "${NEED_TOOLS}" ]] && echo "  tools      :${NEED_TOOLS}"
[[ -n "${NEED_R}" ]]     && echo "  R packages :${NEED_R}"
echo ""

# --- Build the install commands ---------------------------------------------
TOOL_CMD=""
if [[ -n "${NEED_TOOLS}" && -n "${PKG_MGR}" ]]; then
    PKGS=""
    for t in ${NEED_TOOLS}; do
        p="$(pkg_name "${t}" "${PKG_MGR}")"
        case " ${PKGS} " in *" ${p} "*) ;; *) PKGS="${PKGS} ${p}" ;; esac
    done
    case "${PKG_MGR}" in
        brew) TOOL_CMD="brew install${PKGS}" ;;
        apt)  TOOL_CMD="sudo apt-get update && sudo apt-get install -y${PKGS}" ;;
        dnf)  TOOL_CMD="sudo dnf install -y${PKGS}" ;;
        yum)  TOOL_CMD="sudo yum install -y${PKGS}" ;;
    esac
fi

R_CMD=""
if [[ -n "${NEED_R}" ]]; then
    R_LIST="$(echo "${NEED_R}" | tr -s ' ' | sed 's/^ //; s/ $//' | sed "s/ /', '/g")"
    R_CMD="Rscript -e \"install.packages(c('${R_LIST}'), repos='https://cloud.r-project.org')\""
fi

echo "To install:"
[[ -n "${TOOL_CMD}" ]] && echo "    ${TOOL_CMD}"
[[ -n "${R_CMD}" ]]    && echo "    ${R_CMD}"
if [[ -z "${TOOL_CMD}" && -n "${NEED_TOOLS}" ]]; then
    echo "    No supported package manager detected."
    echo "    On an HPC these are usually modules:  module avail | grep -iE 'bcftools|bedtools|htslib|R/'"
fi

if [[ ${HAVE_CONDA} -eq 1 || ${SHOW_CONDA} -eq 1 ]]; then
    echo ""
    echo "Or with conda, which avoids sudo and pins versions per pipeline:"
    echo "    conda env create -f pipelines/methylation/environment.yml"
    echo "    conda env create -f pipelines/sv/environment.yml"
    echo "    conda activate ont-methylation-pipeline   # or ont-sv-pipeline"
fi

if [[ ${DO_INSTALL} -eq 0 ]]; then
    echo ""
    echo "Re-run with --install to have this script run the commands above."
    exit 1
fi

# --- Actually install --------------------------------------------------------
echo ""
echo "--- Installing (--install given) ---"
if [[ "${TOOL_CMD}" == sudo* ]]; then
    echo "  This needs sudo and will prompt for your password."
fi

RC=0
if [[ -n "${TOOL_CMD}" ]]; then
    echo "  \$ ${TOOL_CMD}"
    if ! eval "${TOOL_CMD}"; then
        echo "  [FAIL] tool installation failed" >&2
        RC=1
    fi
fi
if [[ -n "${R_CMD}" ]]; then
    echo "  \$ ${R_CMD}"
    if ! eval "${R_CMD}"; then
        echo "  [FAIL] R package installation failed" >&2
        RC=1
    fi
fi

echo ""
echo "--- Re-checking ---"
# Re-exec rather than duplicating the checks, so there is one source of truth.
exec bash "${BASH_SOURCE[0]:-$0}" ${ONLY_PIPELINE:+--pipeline "${ONLY_PIPELINE}"}
