#!/usr/bin/env bash
# Follow the installation chapter (install.qmd) the way a reader would: run
# every ```bash block for one operating system, in page order, each in a new
# login shell that starts in the home folder (like a new terminal window),
# then check the outputs the chapter says to expect.
#
# Usage: scripts/test-install.sh macos|linux          run the chapter
#        scripts/test-install.sh macos|linux --list   only print the blocks
#
# Running it installs everything into $HOME exactly as the chapter does, so
# only run it on a throwaway machine (CI runner, container, or VM).
set -uo pipefail

# Everything is in main() so that bash reads the whole script before running it.
main() {
  local os=${1:?usage: test-install.sh macos|linux [--list]}
  here=$(cd "$(dirname "$0")" && pwd)
  blocks=$(mktemp -d)
  awk -v os="$os" -v out="$blocks" -f "$here/install-blocks.awk" "$here/../install.qmd" || exit 1

  if [ "${2:-}" = --list ]; then
    for f in "$blocks"/*.sh; do echo "# ---- block $(basename "$f" .sh)"; cat "$f"; done
    exit 0
  fi

  logs=${INSTALL_TEST_LOGS:-$HOME/install-test-logs}
  mkdir -p "$logs"
  case $os in
    macos) shell=(zsh -l -i -c) ;;
    linux) shell=(bash -l -i -c) ;;
  esac

  # Blocks allowed to exit non-zero: the command-line tools are already
  # installed on CI machines, and `gkeyll -v` exits with status 255 even when it
  # works. Their output is still checked below where it matters.
  may_fail='^(xcode-select --install|gkeyll -v)$'

  for f in "$blocks"/*.sh; do
    n=$(basename "$f" .sh)
    first=$(head -n 1 "$f")
    echo "::group::Block $n: $first"
    cat "$f"
    echo "---------------- output ----------------"
    (cd "$HOME" && "${shell[@]}" "set -e; cd ~; $(cat "$f")" < /dev/null) 2>&1 | tee "$logs/$n.log"
    rc=${PIPESTATUS[0]}
    echo "::endgroup::"
    if [ "$rc" -ne 0 ]; then
      if grep -qxE "$may_fail" "$f"; then
        echo "Block $n exited with status $rc, which is expected for it."
      else
        echo "::error title=Install block $n failed::Exit status $rc from block starting: $first"
        exit 1
      fi
    fi
  done

  # Check what the chapter tells readers they should see.
  expect() {  # expect TEXT_IN_BLOCK EXPECTED_OUTPUT
    local f
    f=$(grep -lF -- "$1" "$blocks"/*.sh | head -n 1)
    if [ -z "$f" ]; then
      echo "::error title=Missing block::No block contains '$1'"; return 1
    fi
    if grep -qF -- "$2" "$logs/$(basename "$f" .sh).log"; then
      echo "ok: the block with '$1' printed '$2'"
    else
      echo "::error title=Unexpected output::The block with '$1' did not print '$2'"; return 1
    fi
  }
  status=0
  expect './configure' 'Building APP  moments' || status=1
  expect 'gkeyll -v' 'This is the Gkeyll code' || status=1
  expect 'get_values().shape' '(128, 64, 8)' || status=1
  exit $status
}

main "$@"
