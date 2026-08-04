# =============================================================================
# common/lib_common.sh
# Shared helpers sourced by every pipeline's 00_setup_env.sh.
#
# NOT executable and has no shebang: it is only ever sourced, and it defines
# functions without running anything or setting shell options. The sourcing
# script owns `set`, all output, and all exits.
#
# WHY THIS EXISTS
# The SV and methylation pipelines each grew their own copy of the same config
# reader and tool checker. The copies then diverged, and fixes landed in only
# one of them: the zsh-safe sourcing fix, the '#'-argument guard and the
# duplicate-key guard were all in the methylation copy while the SV pipeline
# still had the original bugs. Adding a third pipeline would have tripled that.
#
# PORTABILITY (see CONTRIBUTING.md)
#   - bash 3.2 compatible: macOS ships 3.2, so no associative arrays,
#     `mapfile`, or ${var,,}.
#   - Safe to source from zsh as well as bash.
#   - No GNU-only tool flags.
# =============================================================================

# -----------------------------------------------------------------------------
# common_config_arg <raw_arg>
# Echoes the config path to use, ignoring an argument that is a shell comment.
#
# zsh does NOT strip inline '#' comments in interactive shells, unlike bash. So
# copying a documented command such as
#     source scripts/bash/00_setup_env.sh   # check the environment
# from a README into a zsh prompt delivers "#" here as $1, and the only symptom
# is a baffling "config file not found: #". No legitimate config path begins
# with '#'.
# -----------------------------------------------------------------------------
common_config_arg() {
    case "${1:-}" in
        '#'*)
            echo "[WARN] Ignoring argument '$1' — that looks like a shell comment." >&2
            echo "       zsh does not strip inline '#' comments interactively." >&2
            echo ""
            ;;
        *) echo "${1:-}" ;;
    esac
}

# -----------------------------------------------------------------------------
# common_yaml_get <config_file> <key>
# First "key: value" line anywhere in the file, quotes and trailing space
# stripped. A targeted grep avoids requiring yq on an HPC.
#
# CONSEQUENCE: returns the FIRST match anywhere, so every key read this way must
# be unique across the file. common_guard_duplicate_keys enforces that rather
# than trusting it.
# -----------------------------------------------------------------------------
common_yaml_get() {
    grep -E "^[[:space:]]*${2}:" "$1" | head -n1 \
        | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//'
}

# -----------------------------------------------------------------------------
# common_guard_duplicate_keys <config_file> <key> [key...]
# Fails if any key appears more than once. A duplicate would be silently
# ignored by the first-match reader above, which is the kind of bug that only
# shows up as a wrong number much later.
# -----------------------------------------------------------------------------
common_guard_duplicate_keys() {
    local cfg="$1"; shift
    local key n dupes=""
    for key in "$@"; do
        n=$(grep -cE "^[[:space:]]*${key}:" "${cfg}")
        [ "${n}" -gt 1 ] && dupes="${dupes} ${key}(x${n})"
    done
    if [ -n "${dupes}" ]; then
        echo "[FATAL] Duplicate config keys in ${cfg}:${dupes}" >&2
        echo "        The config reader returns the first match, so one of these" >&2
        echo "        is being silently ignored. Fix the config before running." >&2
        return 1
    fi
    return 0
}

# -----------------------------------------------------------------------------
# common_resolve_dir <raw_value> <base_dir>
# Echoes an absolute path. Only prepends the base when the value is relative.
#
# Prepending unconditionally leaks test scratch directories into the real repo,
# which is how test logs once ended up in the tracked logs/ directory.
# -----------------------------------------------------------------------------
common_resolve_dir() {
    case "$1" in
        /*) echo "$1" ;;
        *)  echo "$2/$1" ;;
    esac
}

# -----------------------------------------------------------------------------
# common_resolve_filename <config_file> <key> <raw_sample_prefix>
# Reads a filename template and substitutes ${sample.raw_sample_prefix}.
#
# Input filenames use the RAW prefix — the actual Epi2ME output naming — never
# the anonymised SAMPLE_ID, because the raw prefix is what is really on disk.
# -----------------------------------------------------------------------------
common_resolve_filename() {
    common_yaml_get "$1" "$2" | sed "s|\${sample.raw_sample_prefix}|$3|g"
}

# -----------------------------------------------------------------------------
# common_probe_tool <tool>
# Runs `<tool> --version` and reports whether it is usable, not just present.
# Echoes one line to stdout, one of:
#   OK <version-line>
#   FAIL <version-line>    (found on PATH but broken, e.g. missing shared lib)
#   MISS                   (not found on PATH at all)
# Returns 0 for OK, 1 for FAIL/MISS.
#
# PRESENCE IS NOT ENOUGH. A binary can sit on PATH and still be unusable
# through a broken shared library, so the tool is RUN and both its exit status
# and its output are inspected. The SV pipeline once reported [OK] for a
# bcftools that could not load.
#
# This is the single place that logic lives. common_check_tools below and
# scripts/check_dependencies.sh's check_tool() both build on it rather than
# each re-implementing the broken-shared-library detection — that duplication
# is exactly how the SV/methylation copies diverged before common/ existed.
# -----------------------------------------------------------------------------
common_probe_tool() {
    local tool="$1" ver rc
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "MISS"
        return 1
    fi
    ver=$("${tool}" --version 2>&1 | head -n1)
    rc=$?
    case "${ver}" in
        *"cannot open shared object"*|*"error while loading shared libraries"*|*"command not found"*)
            rc=1 ;;
    esac
    if [ "${rc}" -ne 0 ]; then
        echo "FAIL ${ver}"
        return 1
    fi
    echo "OK ${ver}"
    return 0
}

# -----------------------------------------------------------------------------
# common_check_tools <tool_version_log> <tool> [tool...]
# Reports each tool and records its version. Returns 1 if any is missing or
# present-but-broken; the caller decides what to do about it.
# -----------------------------------------------------------------------------
common_check_tools() {
    local log="$1"; shift
    local tool result status ver missing=""

    : > "${log}"
    for tool in "$@"; do
        result="$(common_probe_tool "${tool}")"
        status="${result%% *}"
        ver="${result#* }"
        case "${status}" in
            OK)
                echo "  [OK]   ${tool} -> ${ver}"
                echo "${tool}: ${ver}" >> "${log}"
                ;;
            FAIL)
                echo "  [FAIL] ${tool} found on PATH but failed to run: ${ver}"
                missing="${missing} ${tool}"
                ;;
            MISS)
                echo "  [MISS] ${tool} not found on PATH"
                missing="${missing} ${tool}"
                ;;
        esac
    done

    if [ -n "${missing}" ]; then
        echo ""
        echo "[FATAL] Missing or broken required tools:${missing}"
        echo "        Check what you need with:"
        echo "            bash scripts/check_dependencies.sh"
        echo "        On an HPC these are usually modules:"
        echo "            module avail | grep -iE 'bcftools|bedtools|htslib|R/'"
        return 1
    fi
    return 0
}

# -----------------------------------------------------------------------------
# common_check_writable <dir>
# Verifies dir is actually writable by writing and removing a probe file --
# `-w` reports incorrectly on some network/FUSE mounts. Prints [OK]/[FAIL] to
# stdout/stderr and returns 0/1; the caller decides what to do about it.
#
# Both pipelines' 08_archive_results.sh independently wrote this same
# write-a-probe-file check before common/ existed.
# -----------------------------------------------------------------------------
common_check_writable() {
    local dir="$1" probe
    probe="${dir}/.write_probe_$$"
    if : > "${probe}" 2>/dev/null; then
        rm -f "${probe}"
        echo "  [OK]   ${dir} is writable"
        return 0
    fi
    echo "  [FAIL] ${dir} is not writable" >&2
    return 1
}

# -----------------------------------------------------------------------------
# common_check_free_space <dir> <needed_kb>
# Warns rather than fails if free space cannot be determined (e.g. df flags
# differ on an unusual filesystem) -- an archive should not be blocked by an
# inability to check disk space, only by actually running out of it.
# -----------------------------------------------------------------------------
common_check_free_space() {
    local dir="$1" needed_kb="$2" avail_kb
    avail_kb=$(df -Pk "${dir}" 2>/dev/null | awk 'NR==2 {print $4}')
    if [[ -z "${avail_kb}" ]]; then
        echo "  [WARN] could not determine free space on ${dir}"
        return 0
    fi
    if [[ "${avail_kb}" -lt "${needed_kb}" ]]; then
        echo "  [FAIL] insufficient space: need ~$((needed_kb / 1024)) MB, have $((avail_kb / 1024)) MB" >&2
        return 1
    fi
    echo "  [OK]   space: need ~$((needed_kb / 1024)) MB, have $((avail_kb / 1024)) MB"
    return 0
}

# -----------------------------------------------------------------------------
# common_write_checksums <dest_dir> <checksum_filename> <relative_path> [relative_path...]
# Writes a checksum manifest (recursively, over the given relative paths under
# dest_dir) and immediately verifies it by reading the files back -- a
# checksum file that was only ever generated, never re-read, proves nothing
# about what actually landed on disk.
#
# Prefers sha256sum, falls back to `shasum -a 256`. Echoes one line to stdout:
#   verified:<N>   all N checksums verified
#   failed         generated but did not verify (caller should delete the archive)
#   unavailable    neither sha256sum nor shasum is on PATH
# Returns 0 only for "verified:<N>".
# -----------------------------------------------------------------------------
common_write_checksums() {
    local dest_dir="$1" checksum_name="$2"; shift 2
    local sum_cmd sum_check n

    if command -v sha256sum >/dev/null 2>&1; then
        sum_cmd=(sha256sum); sum_check=(sha256sum -c --quiet)
    elif command -v shasum >/dev/null 2>&1; then
        sum_cmd=(shasum -a 256); sum_check=(shasum -a 256 -c --status)
    else
        echo "unavailable"
        return 1
    fi

    if ! ( cd "${dest_dir}" && find "$@" -type f -print0 | xargs -0 "${sum_cmd[@]}" > "${checksum_name}" ); then
        echo "unavailable"
        return 1
    fi
    n=$(wc -l < "${dest_dir}/${checksum_name}" | tr -d ' ')

    if ( cd "${dest_dir}" && "${sum_check[@]}" "${checksum_name}" ) 2>/dev/null; then
        echo "verified:${n}"
        return 0
    fi
    echo "failed"
    return 1
}

# -----------------------------------------------------------------------------
# common_ref_config <config_file> <repo_dir>
# Resolves reference.config_file to an absolute path.
#
# Read from the pipeline config rather than hardcoded, so a test config can
# point at test references — the SV pipeline originally hardcoded this path and
# the tests could not override it.
# -----------------------------------------------------------------------------
common_ref_config() {
    local raw
    raw="$(common_yaml_get "$1" 'config_file')"
    [ -z "${raw}" ] && { echo ""; return 1; }
    common_resolve_dir "${raw}" "$2"
}
