# These functions allow us to extract the fragment positions directly from the .arrows:
chr_fragments_datatable <- function(h5_con, chr_name) {
    rg <- h5read(h5_con, paste0("/Fragments/",chr_name,"/Ranges"))

    data.table(
        chr = chr_name,
        start = rg[,1],
        end = rg[,1] + rg[,2],
        bc = rep(h5read(h5_con, paste0("/Fragments/",chr_name,"/RGValues")),
                 h5read(h5_con, paste0("/Fragments/",chr_name,"/RGLengths"))),
        n_umi = 1
    )
}

extract_arrow_fragments <- function(arrow_file) {
    h5_con <- H5Fopen(arrow_file)
    contents <- h5ls(h5_con)
    chrs <- contents$name[contents$group == "/Fragments"]

    frags <- chr_fragments_datatable(h5_con, chrs[1])
    for(i in 2:length(chrs)) {
        frags <- rbind(frags, chr_fragments_datatable(h5_con, chrs[i]))
    }

    setorderv(frags, c("chr","start"))

    frags
}

convert_arrow <- function(arrow_file) {
    print(paste("Extracting fragments from",arrow_file))
    fragment_dt <- extract_arrow_fragments(arrow_file)
    fwrite(fragment_dt,
            sub("_archr.arrow","_fragments.tar.gz",arrow_file),
            col.names = FALSE,
            sep = "\t")
}
