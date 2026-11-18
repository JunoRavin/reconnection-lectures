# Extract the ```bash blocks of install.qmd that a reader on one operating
# system would run, in page order, writing one file per block: DIR/01.sh, ...
#
# Usage: awk -v os=macos|linux -v out=DIR -f install-blocks.awk install.qmd
#
# Blocks outside tab sets apply to both systems. Inside a tab set, only the
# block(s) under the "macOS" tab (os=macos) or a tab whose name starts with
# "Linux" (os=linux) are used. Other fenced blocks (sample output) are skipped.

BEGIN {
  if (os != "macos" && os != "linux") {
    print "install-blocks.awk: os must be macos or linux" > "/dev/stderr"
    exit 2
  }
  depth = 0     # nesting depth of ::: divs
  tsdepth = 0   # depth of the enclosing tab set (0 when not in one)
  tab = ""      # current tab name inside that tab set
  incode = 0; keep = 0; n = 0
}

incode {
  if ($0 ~ /^```[ \t]*$/) {
    incode = 0
    if (keep) close(file)
  } else if (keep) {
    print > file
  }
  next
}

/^```/ {
  incode = 1
  keep = 0
  if ($0 ~ /^```bash[ \t]*$/ && (tsdepth == 0 || wanted(tab))) {
    keep = 1
    n++
    file = sprintf("%s/%02d.sh", out, n)
    printf "" > file
  }
  next
}

/^:::+[ \t]*[{A-Za-z]/ {
  depth++
  if ($0 ~ /panel-tabset/) { tsdepth = depth; tab = "" }
  next
}

/^:::+[ \t]*$/ {
  if (depth == tsdepth) { tsdepth = 0; tab = "" }
  if (depth > 0) depth--
  next
}

/^## / && tsdepth > 0 && depth == tsdepth {
  tab = substr($0, 4)
  if (tab != "macOS" && tab !~ /^Linux/)
    print "install-blocks.awk: warning: unknown tab \"" tab "\"; its blocks are skipped" > "/dev/stderr"
  next
}

END {
  if (incode) { print "install-blocks.awk: unterminated code block" > "/dev/stderr"; exit 1 }
  if (n == 0) { print "install-blocks.awk: no bash blocks found" > "/dev/stderr"; exit 1 }
}

function wanted(t) {
  if (os == "macos") return t == "macOS"
  return t ~ /^Linux/
}
