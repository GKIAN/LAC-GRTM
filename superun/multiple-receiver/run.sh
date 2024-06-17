#!/bin/bash

#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Mon 17 Jun 2024 04:53:06 PM CST
#-------------------------------------------------------------------------------

set -e

rcvin=./receivers.in
path=../../GRTM-E3D

export OMP_NUM_THREADS=4

exe=${path}/bin/grtcsgram
demoinput=${path}/demo/input.conf
rowpref="coordinate_receiver"

n=`cat ${rcvin} | wc -l`
for ((i = 0; i < ${n}; i++))
do
  rxyz=`sed -n "$((${i} + 1))p" ${rcvin}`
  cat ${demoinput} | sed "s/^${rowpref} =.*/${rowpref} = ${rxyz}/" > input.conf
  sed -i "s/^output_prefix =.*/output_prefix = 'out-${i}'/" input.conf

  echo ">> Running [${exe}] for the receiver #${i} ..."
  ${exe} input.conf
done
