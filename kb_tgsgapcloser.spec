/*
A KBase module: kb_tgsgapcloser
Fill N-gaps in scaffolded genome assemblies with long reads (or pre-assembled
contigs) using TGS-GapCloser v1.2.1, with optional Racon correction of the
candidate fill sequences.
*/

module kb_tgsgapcloser {

    /* A workspace object reference, e.g. 123/4/5 */
    typedef string obj_ref;

    /*
    workspace_name - workspace to save the gap-filled Assembly and report in
    assembly_ref - scaffolded Assembly containing N-gaps
    reads_ref - SingleEndLibrary of long reads, or an Assembly of contigs
    output_assembly_name - name of the gap-filled Assembly
    read_type - ont | pb_clr | pb_hifi
    error_correction - racon | none
    racon_rounds - Racon iterations on candidate fills (1-4)
    minimap2_preset - auto | asm5 | asm10 | asm20
    min_idy, min_match - candidate filters; blank means the tool's
        read-type default (ont 0.3 / 300, pb 0.2 / 200)
    min_nread - minimum reads that must bridge a gap
    max_candidate - maximum candidate fills per gap
    chunk - number of chunks for candidate correction
    gapsize_check - 1 to reject fills that disagree with the gap size
    threads - blank means all available CPUs
    */
    typedef structure {
        string workspace_name;
        obj_ref assembly_ref;
        obj_ref reads_ref;
        string output_assembly_name;
        string read_type;
        string error_correction;
        int racon_rounds;
        string minimap2_preset;
        float min_idy;
        int min_match;
        int min_nread;
        int max_candidate;
        int chunk;
        int gapsize_check;
        int threads;
    } TGSGapCloserParams;

    typedef structure {
        string report_name;
        string report_ref;
        obj_ref assembly_ref;
        UnspecifiedObject gapfill_summary;
    } TGSGapCloserOutput;

    funcdef run_tgsgapcloser(TGSGapCloserParams params)
        returns (TGSGapCloserOutput output) authentication required;
};
