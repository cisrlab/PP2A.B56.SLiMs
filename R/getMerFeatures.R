
getMerFeatures<-function(Sequences, min_mer = 4, max_mer=12, filter.sd=TRUE, filter.redundant=TRUE) {

    seq_counts = zcurve$getFeatureCounts(Sequences, min_mer=min_mer, max_mer = max_mer)

    seq_fs = matrix(0, nrow=length(Sequences), ncol=nrow(seq_counts));
    colnames(seq_fs) = seq_counts$seq;
    for (seq_idx in 1:nrow(seq_fs)) {
        seq_counts = zcurve$getFeatureCounts(Sequences[seq_idx], min_mer=min_mer, max_mer=max_mer)
        seq_fs[seq_idx,seq_counts$seq] = seq_counts$count;
    }
    if (filter.sd) {
        sds =matrixStats::colSds(seq_fs);
        seq_fs = seq_fs[,sds != 0];
    }
    seq_fs = as.data.frame(seq_fs, stringsAsFactors=FALSE)

    seq_fn = colnames(seq_fs);

    keep = rep(TRUE, length(seq_fn))
    if (filter.redundant) {
        for (idx1 in 1:(length(seq_fn)-1)) {
            if (keep[idx1]) {
                v1_vec = seq_fs[,idx1];
                for (idx2 in (idx1+1):length(seq_fn)) {
                    v2_vec = seq_fs[,idx2];
                    if (sum(v1_vec == v2_vec) == nrow(seq_fs)) {
                        keep[idx2] = FALSE;
                    }
                }
            }
        }
    }

    seq_fs = seq_fs[,keep, drop=FALSE];
    #rownames(seq_fs) = Sequences;
    return(seq_fs);

}
