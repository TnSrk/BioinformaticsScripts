more | awk '{if ($0 ~ /^>/) print "#"$0"#";else print $0"#\t"length($0) }' | tr -d '\n' | sed "s/#>/\n>/g" | awk '{if ($2 < 200 ) print $0}' | awk '{print $1}' | tr '#' '\n' | grep -v "^$"
