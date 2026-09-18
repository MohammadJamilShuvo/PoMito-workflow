#!/usr/bin/env bash
set -euo pipefail
REFROOT="${1:-${POMITO_MITOS_REFDIR:-${HOME}/.pomito/mitos_refdata}}"
REFSEQVER="${POMITO_MITOS_REFSEQVER:-refseq89m}"
ARCHIVE="${REFSEQVER}.tar.bz2"
URL="https://zenodo.org/records/4284483/files/${ARCHIVE}?download=1"
EXPECTED_MD5="bb6325e27612e61a6995e88e8ffcecf8"
mkdir -p "${REFROOT}"
if [[ -d "${REFROOT}/${REFSEQVER}" && -f "${REFROOT}/${REFSEQVER}/auxinfo.json" ]]; then
  echo "[OK] MITOS2 reference data already present: ${REFROOT}/${REFSEQVER}"; exit 0
fi
tmp="${REFROOT}/${ARCHIVE}"
echo "[INFO] Downloading MITOS2 ${REFSEQVER} reference data"
if command -v curl >/dev/null 2>&1; then curl -L --fail --retry 3 -o "${tmp}" "${URL}"
elif command -v wget >/dev/null 2>&1; then wget -O "${tmp}" "${URL}"
else echo "[ERROR] curl or wget is required." >&2; exit 1; fi
actual_md5="$(md5sum "${tmp}" | awk '{print $1}')"
[[ "${actual_md5}" == "${EXPECTED_MD5}" ]] || { echo "[ERROR] MITOS2 reference-data MD5 mismatch: ${actual_md5}" >&2; exit 1; }
tar -xjf "${tmp}" -C "${REFROOT}"; rm -f "${tmp}"
[[ -d "${REFROOT}/${REFSEQVER}" ]] || { echo "[ERROR] Expected directory not created: ${REFROOT}/${REFSEQVER}" >&2; exit 1; }
echo "[OK] MITOS2 reference data ready: ${REFROOT}/${REFSEQVER}"
