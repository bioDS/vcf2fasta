library("beter")

if(sys.nframe() == 0){
    args = commandArgs(TRUE)
    beter::process_template(
        template=args[1], output=args[3], alignment=args[2],
        parameters=list("datatype"="nucleotideDiploid16"))
    }
