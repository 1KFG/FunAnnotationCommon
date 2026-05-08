#!/usr/bin/env bash
# Delete annotate/<folder> directories where predict_misc contains only ab_initio_parameters,
# indicating a failed/incomplete funannotate predict run with nothing worth keeping.
#
# Usage:
#   cleanup_incomplete_predict.sh [OPTIONS] [FOLDER...]
#   cat folder_list.txt | cleanup_incomplete_predict.sh [OPTIONS]
#
# Options:
#   -a, --annotate DIR   Path to annotate/ directory (default: annotate/ relative to script)
#   -n, --dry-run        Print what would be deleted without deleting anything
#   -d, --debug          Verbose output for every folder checked
#   -h, --help           Show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANNOTATE_DIR="$(realpath "annotate")"
DRY_RUN=0
DEBUG=0

usage() {
    grep '^#' "$0" | sed 's/^# \?//' | tail -n +2
    exit 0
}

log_debug() { [[ $DEBUG -eq 1 ]] && echo "[DEBUG] $*" >&2 || true; }
log_info()  { echo "[INFO]  $*"; }
log_warn()  { echo "[WARN]  $*" >&2; }

# ── argument parsing ──────────────────────────────────────────────────────────
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -a|--annotate)  ANNOTATE_DIR="$(realpath "$2")"; shift 2 ;;
        -n|--dry-run)   DRY_RUN=1; shift ;;
        -d|--debug)     DEBUG=1; shift ;;
        -h|--help)      usage ;;
        --) shift; POSITIONAL+=("$@"); break ;;
        -*) echo "Unknown option: $1" >&2; exit 1 ;;
        *)  POSITIONAL+=("$1"); shift ;;
    esac
done

if [[ ! -d "$ANNOTATE_DIR" ]]; then
    echo "ERROR: annotate directory not found: $ANNOTATE_DIR" >&2
    exit 1
fi

log_debug "annotate dir : $ANNOTATE_DIR"
log_debug "dry-run      : $DRY_RUN"

[[ $DRY_RUN -eq 1 ]] && log_info "DRY-RUN mode — nothing will be deleted"

# ── build list of folders to check ───────────────────────────────────────────
# Prefer positional args; fall back to stdin if none given and stdin is a pipe.
get_folders() {
    if [[ ${#POSITIONAL[@]} -gt 0 ]]; then
        printf '%s\n' "${POSITIONAL[@]}"
    elif [[ ! -t 0 ]]; then
        cat
    else
        echo "ERROR: provide folder names as arguments or via stdin" >&2
        exit 1
    fi
}

# ── main loop ─────────────────────────────────────────────────────────────────
deleted=0
skipped=0
missing=0

get_folders $ANNOTATE_DIR

while IFS= read -r folder; do
    [[ -z "$folder" ]] && continue
    log_debug "in loop: -> $folder"

    target="${ANNOTATE_DIR}/${folder}"
    misc="${target}/predict_misc"

    if [[ ! -d "$target" ]]; then
        log_warn "folder not found, skipping: $target"
        (( missing++ )) || true
        continue
    fi

    if [[ ! -d "$misc" ]]; then
        log_debug "$folder: no predict_misc — skipping"
        (( skipped++ )) || true
        continue
    fi

    # Count entries in predict_misc
    mapfile -t entries < <(ls -A "$misc")
    count=${#entries[@]}

    log_debug "$folder: predict_misc has $count item(s): ${entries[*]:-<empty>}"

    if [[ $count -eq 1 && "${entries[0]}" == "ab_initio_parameters" && -d "${misc}/ab_initio_parameters" ]]; then
        if [[ $DRY_RUN -eq 1 ]]; then
            log_info "WOULD DELETE: $target"
        else
            log_info "Deleting: $target"
            rm -rf "$target"
        fi
        (( deleted++ )) || true
    else
        log_debug "$folder: predict_misc contents don't match criteria — skipping"
        (( skipped++ )) || true
    fi
done < <(get_folders)

echo ""
echo "── Summary ──────────────────────────────────────"
[[ $DRY_RUN -eq 1 ]] && echo "  Would delete : $deleted"   \
                      || echo "  Deleted      : $deleted"
echo "  Skipped      : $skipped"
echo "  Not found    : $missing"
