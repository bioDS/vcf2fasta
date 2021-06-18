library("phyloRNA")
library("beter")
library("argparser")


parse_args = function(){
    parser = argparser::arg_parser(
        "Calculate a few summary statistics for a fasta file"
        )
    parser = add_argument(parser, "fasta", type="character",
        help="An input fasta file")
    parser = add_argument(parser, "stats", type="character",
        help="An output path for a file with calculated statistics")
    parser = add_argument(parser, "--filtered", type="character",
        help="Path for a filtered fasta file with constant sites removed")
    args = argparser::parse_args(parser)
    subset(args, !names(args) %in% c("", "help", "opts"))
    }


main = function(fasta, stats, filtered=NULL){
    seq = beter::read_fasta(fasta)
    tab = phyloRNA::seq2tab(seq)
    tab_nc = phyloRNA::remove_constant(tab, margin=2, unknown="?")

    text = paste0(
        "Taxa: ", nrow(tab), "\n",
        "Sites: ", ncol(tab), "\n",
        "Constant sites: ", ncol(tab) - ncol(tab_nc), "\n",
        "Unknown sites (%): ",
        round( sum(tab == "?")/prod(dim(tab) ), 2), "\n",
        "Unknown sites (%) (no constant): ",
        round( sum(tab_nc == "?")/prod(dim(tab_nc)), 2), "\n"
        )
    writeLines(text, stats)
    
    if(!is.null(filtered)){
        seq_nc = phyloRNA::tab2seq(tab_nc)
        beter::write_fasta(seq_nc, filtered)        
        }
    }


if(sys.nframe() == 0){
    args = parse_args()
    do.call(main, args)
    }
