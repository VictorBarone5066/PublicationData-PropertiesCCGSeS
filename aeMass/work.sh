#!/bin/bash

phases=("Tetra" "Ortho")
concs=("0000" "0125" "0250" "0375" "0500" "0625" "0750" "0875" "1000")

sx="-x"
aem=".aem"

for ph in "${phases[@]}"; do
	for co in "${concs[@]}"; do
		nameEig="eig$ph$sx$co$aem"
		nameIn="in$ph$sx$co$aem"
		nameOut="out$ph$sx$co$aem"
		echo "$nameIn $nameEig $nameOut"

		./aem -i "$nameIn" -e "$nameEig" -o "$nameOut"
	done
done

