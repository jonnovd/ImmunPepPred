include { 
    PREP_DEEPTAP_INPUT 
    RUN_DEEPTAP
} from '../processes/tap_processes'

workflow TAP_WORKFLOW {
    main:
        dtin_ch = channel.fromPath(params.deeptap_in)
        PREP_DEEPTAP_INPUT(dtin_ch)

        PREP_DEEPTAP_INPUT.out.view()

        RUN_DEEPTAP(PREP_DEEPTAP_INPUT.out)

    emit:
        out = RUN_DEEPTAP.out
}