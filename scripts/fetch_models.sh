#!/bin/bash
# Fetch MuJoCo Menagerie robot models for cross-validation
#
# Downloads MJCF files for well-known robots into models/menagerie/
# using sparse checkout to avoid cloning the entire repository.
#
# Usage: ./scripts/fetch_models.sh

set -euo pipefail

MODELS_DIR="$(cd "$(dirname "$0")/.." && pwd)/models/menagerie"
TMP_REPO="/tmp/mujoco_menagerie_$$"

# Robots to download (from google-deepmind/mujoco_menagerie)
ROBOTS=(
    "unitree_go1"
    "unitree_go2"
    "unitree_h1"
    "unitree_g1"
)

mkdir -p "$MODELS_DIR"

echo "=== Cloning MuJoCo Menagerie (sparse checkout) ==="
git clone --depth 1 --filter=blob:none --sparse \
    https://github.com/google-deepmind/mujoco_menagerie.git "$TMP_REPO" 2>&1

cd "$TMP_REPO"
git sparse-checkout set "${ROBOTS[@]}" 2>&1

for robot in "${ROBOTS[@]}"; do
    if [ -d "$robot" ]; then
        echo "  Copying $robot..."
        cp -r "$robot" "$MODELS_DIR/"
        # Remove .git directory if accidentally included
        rm -rf "$MODELS_DIR/$robot/.git" 2>/dev/null || true
    else
        echo "  WARNING: $robot not found in repository"
    fi
done

cd /
rm -rf "$TMP_REPO"

echo ""
echo "=== Models downloaded to $MODELS_DIR ==="
ls -la "$MODELS_DIR"/*/ 2>/dev/null || echo "  (empty)"
echo ""
echo "Done."
