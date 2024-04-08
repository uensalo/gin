#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 file1.ging file2.ging"
  exit 1
fi

awk -v file1="$1" -v file2="$2" '
  BEGIN {
    FS=OFS="\t"
  }
  FNR == NR { # Processing the first file
    if (/^V/ || /^E/) unique[$0] = 1
    next
  }
  { # Processing the second file
    if ((/^V/ || /^E/) && $0 in unique) delete unique[$0]
    else if (/^V/ || /^E/) diff[$0] = 1
  }
  END {
    for (line in unique) {
      if (line ~ /^V/) {
        split(line, arr, FS)
        print "V", arr[2], "in", file1, "is not in", file2
      } else if (line ~ /^E/) {
        split(line, arr, FS)
        print "E", arr[2], " ", arr[3], "in", file1, "is not in", file2
      }
    }
    for (line in diff) {
      if (line ~ /^V/) {
        split(line, arr, FS)
        print "V", arr[2], "in", file2, "not in", file1
      } else if (line ~ /^E/) {
        split(line, arr, FS)
        print "E", arr[2], " ", arr[3], "in", file2, "is not in", file1
      }
    }
  }
' "$1" "$2"