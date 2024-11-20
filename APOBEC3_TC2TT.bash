#grep -e "      .       C       T       " -e "  .       T       <NON_REF>       "  gatk4_OUT.filtered.normalized.snpeffed.vcf \
grep -e "       .       C       T       " -e "  .       T       <NON_REF>       "  $1 \
        | grep -A 1 "   .       T       <NON_REF>       " | grep -B 1 " .       C       T       "  | grep -A 1 "        .       T       <NON_REF>       "\
        | awk -v chk=0 -v cmp=0 -v tmps="" '{if (chk==0 && cmp==0) {cmp=$2;chk=1;tmps=$0} else { if (cmp != 0 && $2 - cmp == 1) {print $1} }  ;if ($1 == "--") {chk=0;cmp=0} }'\
        | awk '{match($10, /,[0-9]+:[0-9]+:/); $10=substr($10, RSTART + 1, RLENGTH - 2);print  $10 "\t" $0}' \
        | sed "s/:/\t/" \
        | awk '{if ($2 > 0) {print $1/$2 "\t" $3 ":" $4}}'
