############### renum
cd /mnt/c/Users/lzhang/OneDrive\ -\ AbacusBio Ltd/Project_internal/Potato\ mash/blupf90/vce/renum_gbs
time renumf90.exe renum.par |tee renum.log
# then manually create a XrefID for haplotype and for gbs a22
# Run make_xrefid.R

############### obtain Gi matrices
cd ../pregs_hap # "gi"
time preGSf90.exe pre.par |tee pre.log # The Gi is cleaned
cd ../pregs_gbs_a22 # "a22i"
time preGSf90.exe pre.par |tee pre.log # The Gi is cleaned
cd ../pregs_gbs # "ai"
time preGSf90.exe pre.par |tee pre.log # The Gi is cleaned
# cp Gi Gi_pre
# time preGSf90.exe pre_step1.par |tee pre_1.log # clean
# time preGSf90.exe pre_step2.par |tee pre_2.log # save Gi/.

# then make H inverse from make_Hinv.R

############## vce and ge of full data
cd ../vce_gbs
time blupf90+.exe ssgblup.par |tee vce.log
cd ../vce_hap
time blupf90+.exe ssgblup.par |tee vce.log # check if VCE different
cd ../vce_blend
time blupf90+.exes sgblup.par --dense |tee vce.log # check if VCE different

############## ge
cd renum_gbs
# remove year 17 data. Run remove_young_record.R
cd ../ge_gbs_reduced
time blupf90+.exe ssgblup.par |te e ge.log
