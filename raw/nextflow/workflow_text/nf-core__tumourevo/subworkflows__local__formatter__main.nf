//
// FORMATTING SUB-WORKFLOW
//

include { CNAQC2TSV } from '../../../modules/local/cnaqc2tsv/main'
include { CNA2CNAQC } from '../../../modules/local/cna2cnaqc/main'
include { VCF2CNAQC } from '../../../modules/local/vcf2cnaqc/main'


workflow FORMATTER {
    take:
        input
        extension

    main:

        if (extension == "vcf"){
                VCF2CNAQC(input)
                out = VCF2CNAQC.out.rds

        } else if (extension == "cna"){
                input_cna = input.map{ meta, cna  ->
                    [ meta, cna[0], cna[1] ]
                }
                CNA2CNAQC(input_cna)
                out = CNA2CNAQC.out.rds

        } else if (extension == "rds"){ // for pyclone-vi
                CNAQC2TSV(input)
                out = CNAQC2TSV.out.tsv
        }

    emit:
        out
}
