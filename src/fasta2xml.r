library("beter")

if(sys.nframe() == 0){
    args = commandArgs(TRUE)

    alignment_id = basename(tools::file_path_sans_ext(args[3]))
    beter::process_template(
        template = args[1],
        output = args[3],
        alignment = args[2],
        parameters = list(
            "datatype" = "nucleotideDiploid16",
            "alignment_id" = alignment_id
            )
        )
    }
