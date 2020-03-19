# ANOMALY Pipeline bash mode

cd ~/Repository/LRF/00_erifa_bak/ranomaly_test

# DADA2
Rscript ~/Repository/LRF/ranomaly/scripts/dada2_process.R -p reads -c TRUE -q TRUE

# Assign
# assign_taxo(dada_res = dada_res, id_db = "/home/erifa/bank/silva/SILVA_SSU_r132_March2018.RData")

Rscript ~/Repository/LRF/ranomaly/scripts/assign_taxo.R -r ./dada2_out/robjects.Rdata -i /home/erifa/bank/silva/SILVA_SSU_r132_March2018.RData
