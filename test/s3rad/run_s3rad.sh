#!/usr/bin/env bash
set -u

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)
# S3RAD_DECK lets a pilot run use a DERIVED deck (see mkpilotdeck.py).
# The committed deck's hash is then not the right thing to check, so the
# check is replaced by recording the derived deck's own hash in
# provenance.txt together with the fact that it is not the committed one.
DECK=${S3RAD_DECK:-"$SCRIPT_DIR/m12_s3rad_gc24_w.inp"}
EXE=${CCX_EXE:-"$ROOT/src/ccx_2.23"}
RUN_DIR=${1:-"$SCRIPT_DIR/_runs/$(date +%Y%m%d-%H%M%S)"}
if [ "$#" -gt 0 ]; then shift; fi

case "$EXE" in
    /*) ;;
    *) EXE=$(CDPATH= cd -- "$(dirname -- "$EXE")" && pwd)/$(basename -- "$EXE") ;;
esac
case "$RUN_DIR" in
    /*) ;;
    *) RUN_DIR="$PWD/$RUN_DIR" ;;
esac

EXPECTED_DECK_SHA=2fb0cf4e3554282e1f85cf641c788939bc7a2fa5918dd842f54d38844dfa5391

if [ ! -f "$DECK" ]; then
    echo "s3rad deck not found: $DECK" >&2
    exit 2
fi
if [ ! -x "$EXE" ]; then
    echo "CalculiX executable is not executable: $EXE" >&2
    echo "Set CCX_EXE to a PARDISO-enabled ccx_2.23 binary." >&2
    exit 2
fi

ACTUAL_DECK_SHA=$(sha256sum "$DECK" | awk '{print $1}')
DECK_IS_COMMITTED=yes
if [ "$ACTUAL_DECK_SHA" != "$EXPECTED_DECK_SHA" ]; then
    if [ -n "${S3RAD_DECK:-}" ]; then
        DECK_IS_COMMITTED=no
        echo "NOTE: running a DERIVED deck, not the committed one." >&2
        echo "      $DECK" >&2
        echo "      sha256 $ACTUAL_DECK_SHA" >&2
    else
        echo "deck hash mismatch" >&2
        echo "expected: $EXPECTED_DECK_SHA" >&2
        echo "actual:   $ACTUAL_DECK_SHA" >&2
        exit 2
    fi
fi

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-6}
export MKL_CBWR=${MKL_CBWR:-COMPATIBLE}
export CCX_DAMAGE_AUTOSPC=1.e-3
export CCX_DAMAGE_DEADALL=1.e-2
export CCX_DAMAGE_DELETE_MAT=ALL
export CCX_DAMAGE_LINESEARCH=ADAPTIVE
export CCX_DAMAGE_REEQ_RESCUE2=1
export CCX_DAMAGE_REEQ_SCALE=PHYSICAL
export CCX_DAMAGE_TANGENT=UNSYM
export CCX_DAMAGE_TOPOLOGY=DEFERRED
export CCX_DAMAGE_TR_DOGLEG=1
export CCX_DAMAGE_VISCOSITY=1.e-4
export CCX_FRACTURE_TERMINATION=FACE_X0_NSET:FACE_XL_NSET
export CCX_PARDISO_REUSE_SYMBOLIC=1

unset CCX_DISSIPATION_CONTROL
unset CCX_DISSIPATION_TARGET
unset CCX_FRACTURE_LINK
# CCX_FRACTURE_DEADFACET is retired.  Unsetting it here is what kept this
# runner - the only one that armed the connectivity test - from ever seeing
# the specimen come apart, so it is worth saying plainly rather than deleting
# the line silently: src/loadpath.c now applies the judgement unconditionally.
unset CCX_FRACTURE_DEADFACET

for kv in "$@"; do
    name=${kv%%=*}
    value=${kv#*=}
    case "$name" in
        ''|*[!A-Za-z0-9_]*)
            echo "invalid environment override: $kv" >&2
            exit 2
            ;;
    esac
    if [ "$value" = "$kv" ]; then
        echo "override must be NAME=VALUE: $kv" >&2
        exit 2
    elif [ -z "$value" ]; then
        unset "$name"
    else
        export "$kv"
    fi
done

mkdir -p "$RUN_DIR"
cp "$DECK" "$RUN_DIR/m.inp"

EXE_SHA=$(sha256sum "$EXE" | awk '{print $1}')
{
    echo "deck=$DECK"
    echo "deck_sha256=$ACTUAL_DECK_SHA"
    echo "deck_is_committed=$DECK_IS_COMMITTED"
    echo "executable=$EXE"
    echo "executable_sha256=$EXE_SHA"
    echo "started=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "overrides=$*"
    echo "--- environment ---"
    env | grep -E '^(CCX_|OMP_|MKL_)' | sort
} > "$RUN_DIR/provenance.txt"

echo "Run directory: $RUN_DIR"
echo "Deck SHA-256: $ACTUAL_DECK_SHA"
echo "Executable SHA-256: $EXE_SHA"
echo "The committed deck requests PARDISO."

start=$(date +%s)
(
    cd "$RUN_DIR" || exit 2
    "$EXE" -i m > run.log 2>&1
)
rc=$?
finish=$(date +%s)

{
    echo "return_code=$rc"
    echo "wall_seconds=$((finish-start))"
} >> "$RUN_DIR/provenance.txt"

echo "Return code: $rc"
echo "Wall seconds: $((finish-start))"
tail -n 5 "$RUN_DIR/m.sta" 2>/dev/null || true
exit "$rc"
