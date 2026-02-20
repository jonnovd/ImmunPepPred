include { 
    PREP_DEEPTAP_INPUT 
    RUN_DEEPTAP
} from '../processes/tap_processes'

workflow TAP_WORKFLOW {
    take:
        deeptap_in
    
    main:
        //dtin_ch = channel.fromPath(deeptap_in)
        //PREP_DEEPTAP_INPUT(dtin_ch)
        PREP_DEEPTAP_INPUT(deeptap_in)

        RUN_DEEPTAP(PREP_DEEPTAP_INPUT.out)

    emit:
        out = RUN_DEEPTAP.out
}