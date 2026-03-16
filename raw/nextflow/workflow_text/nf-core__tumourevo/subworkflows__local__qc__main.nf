//
// QC SUB-WORKFLOW
//

include { TINC } from '../../../modules/nf-core/tinc/main'
include { CNAQC } from '../../../modules/nf-core/cnaqc/main'
include { JOIN_CNAQC } from '../../../modules/local/join_cnaqc/main'


workflow QC {
    take:
        input

    main:
        TINC(input)

        input_cnaqc = input.map{meta, cna, snv ->
                def sample =  meta.tumour_sample
                [meta, snv, cna, sample]
        }

        CNAQC(input_cnaqc)
        in_join_cnaqc = CNAQC.out.qc_rds.map{ meta, rds ->
            def sample = meta.tumour_sample
            meta = meta + [id: "${meta.dataset}_${meta.patient}"]
            [meta.subMap('dataset', 'patient', 'id'), rds, sample]}
            .groupTuple()

        JOIN_CNAQC(in_join_cnaqc)

    emit:
        rds_cnaqc = CNAQC.out.qc_rds
        plot_cnaqc_data_rds = CNAQC.out.data_plot_rds
        plot_cnaqc_qc_rds = CNAQC.out.qc_plot_rds
        plot_cnaqc_data = CNAQC.out.plot_pdf_data
        plot_cnaqc_qc = CNAQC.out.plot_pdf_qc

        plot_rds_tinc = TINC.out.plot_rds
        rds_tinc = TINC.out.rds
        pdf_tinc = TINC.out.plot_pdf
        csv_tinc = TINC.out.tinc_csv

        join_cnaqc_ALL = JOIN_CNAQC.out.rds_all
        join_cnaqc_PASS = JOIN_CNAQC.out.rds_pass
}
