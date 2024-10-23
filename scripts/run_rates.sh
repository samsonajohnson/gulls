#!/bin/bash

grmatch -r gulls_surot2d_H2023.sources -i gulls_surot2d_H2023.lenses -o - --match-coords --col-ref 2,3 --col-inp 2,3 | while read sidx sl sb sdl sdb ssa sfile lidx ll lb ldl ldb lsa lfile; do
	echo -n "$sidx $sl $sb $sdl $sdb $ssa $sfile $lidx $ll $lb $ldl $ldb $lsa $lfile "
	python ~/gulls_fz/scripts/calc_rate.py $sfile $ssa $lfile $lsa
done
