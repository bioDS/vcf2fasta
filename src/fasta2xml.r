library("beter")
library("argparser")


args_parser = function(){
    parser = arg_parser("Produce a BEAST XML from a beter XML template and a fasta file.")
    parser = add_argument(parser, "template", type="character", help="A `beter` XML template")
    parser = add_argument(parser, "fasta", type="character", help="A fasta file")
    parser = add_argument(parser, "output", type="character", help="Beast XML file with sequences")
    parser = add_argument(parser, "--datatype", default="standard",
        help="Beast sequence data type.")
    args = parse_args(parser)
    subset(args, !names(args) %in% c("", "help", "opts"))
    }

main = function(template, fasta, output, datatype="standard"){
    alignment_id = basename(tools::file_path_sans_ext(output))
    beter::process_template(
        template, output, alignment=fasta,
        parameters = list("datatype" = datatype, "alignment_id" = alignment_id)
        )
    }

if(sys.nframe() == 0){
    args = args_parser()
    do.call(main, args)
    }
