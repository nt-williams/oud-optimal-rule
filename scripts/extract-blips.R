# load functions
box::use(../R/blip, here[here], glue[glue])

blip_type <- "type1"

# importing estimated TSMs
# tsms <- readRDS(here("data", "drv", "onestep-tsm-imputed.rds"))
tsms <- readRDS(here("data", "drv", "onestep-tsm-imputed-no27bup.rds"))

# extract blip
cates <- lapply(tsms, function(x) {
  switch(blip_type,
         type1 = blip$type_1_blip(x$met, x$nal, ref = x$bup),
         type2 = blip$type_2_blip(x$met, x$bup, x$nal))
})

saveRDS(cates, here("data", "drv", glue("{blip_type}-blips-no27bup.rds")))
