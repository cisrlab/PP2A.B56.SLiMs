
zcurve=list();

zcurve$dna.base = c("A","G","C","T");
zcurve$aa.base = c("G", "A", "V", "L", "M", "I", "F", "Y", "W", "K", "H", "N", "S", "T", "C", "P", "D", "Q", "E");




zcurve$getFeatureBase<-function(mer_size=4, base = zcurve$dna.base) {

  if (mer_size == 1) {
     return(base);
  } else {
    ans = c();
    for (c in base) {
      ans = c(ans, paste0(c, zcurve$getFeatureBase(mer_size-1, base)));
    }
    return(ans);
  }
}


zcurve$getFeatureLabels<-function(mer_size=4, base = zcurve$dna.base) {

  ans = c();

  for (m in 1:mer_size) {
    ans = c(ans,zcurve$getFeatureBase(m, base));
  }

  return(unique(ans));

}

zcurve$getFeatureCounts<-function(sequences, min_mer=1, max_mer=4, debug=FALSE) {

  #mer_size = min(nchar(sequence), mer_size);
  
#  f.names = zcurve$getFeatureLabels(mer_size, base = base);
 
#  ans = rep(0,length(f.names));
#  names(ans) = f.names;
  ans = c();
  seqs = unique(sequences);
  for (sequence_idx in 1:length(seqs)) {
    sequence = seqs[sequence_idx];
    #if (debug){cat("idx:",sequence_idx," ",length(seqs)," ",sequence,"\n");}
    for (m in min_mer:max_mer) {
      if ((nchar(sequence)-m+1) >= 1) {
        for (i in 1:(nchar(sequence)-m+1)) {
          s = substr(sequence, i, i+m-1);
          #print(s);
          ans = c(ans,s);
        }
      }
    }
  }
  ans.tbl = table(ans);
  
  ans.df = as.data.frame(ans.tbl);
  colnames(ans.df) = c("seq","count")
  ans.df$seq = as.character(ans.df$seq);
  ans.df$length = nchar(ans.df$seq);
  return(ans.df);
}

zcurve$bindCounts<-function(counts1, counts2) {
  ans = merge(counts1, counts2, by.x="seq",by.y="seq", all.x=TRUE, all.y=TRUE);
  NAs = is.na(ans);
  ans[NAs] = 0;
  return(ans);
        
}

zcurve$createDF<-function(seqs, mer_size=4) {
  useqs = unique(seqs);
  ans.df = zcurve$getFeatureCounts(useqs[1], mer_size=mer_size);
  if (length(useqs) > 1) {
    for (seq.idx in 2:length(useqs)) {
      seq = useqs[seq.idx];
      ans.df = zcurve$bindCounts(ans.df, zcurve$getFeatureCounts(seq, mer_size=mer_size));
    }
  }
  
  colnames(ans.df)[2:(length(useqs)+1)] = useqs;

  rownames(ans.df) = ans.df[,1];
  ans.df = ans.df[,-1];

  return(ans.df);

}

zcurve$getPossibleMers<-function(sequences, mer.size=4, filter=TRUE) {
  cat("Finding bases\n");
  bases = unique(unlist(strsplit(sequences,""))) #Find all unique letters in the list of sequences
  print(bases);
  cat("Generating\n");
  mers = gtools::permutations(n=length(bases),r=mer.size,v=bases,repeats.allowed=T) #Generate all possible kmer sequences using the letters as base
  
  cat("Concatinating ",nrow(mers),"\n");
  mers = apply(mers, 1, function(l){paste0(l,collapse="")})
  
  if (filter) {
    cat("Filtering ",length(mers),"\n");
    library(operators)
    keep = apply(as.matrix(mers), 1, function(l) {sequences %~+% l})
    mers = mers[keep];
  }
  return(mers);
}

zcurve$getBases<-function(sequences) {
  bases = unique(unlist(strsplit(sequences,""))) #Find all unique letters in the list of sequences
  return(bases);
}

zcurve$findMersR<-function(sequences, prev_mer="", bases = zcurve$getBases(sequences), max.mer=4) {
  ans = c();
  if (length(sequences) == 0 || nchar(prev_mer) == max.mer) {
    return(ans);
  }
  for (base in bases) {
    new_mer = paste0(prev_mer,base);
    new_sequences = sequences[grep(new_mer, sequences)]
    if (length(new_sequences) > 0) {
      ans = c(ans, new_mer);
      ans = c(ans, zcurve$findMersR(new_sequences,new_mer,bases, max.mer=max.mer))
    }
  }
  return(ans);
  
}

zcurve$findMers<-function(sequences, bases = zcurve$getBases(sequences), min.mer = 4, max.mer = 4) {
  ans = zcurve$findMersR(sequences, bases=bases, max.mer=max.mer);
  ans = ans[nchar(ans) >= min.mer];
  return(ans);
}

zcurve$writeMers<-function(sequences, prev_mer="", bases=zcurve$getBases(sequences), min.mer = 4, max.mer=4, file.cons = list()) {
  #cat("writeMers:start()", " prev_mer:",prev_mer,"\n");
  #print(names(file.cons));
  
  if (prev_mer == "") { # If 1st call, make sure sequences are unique.
    sequences = unique(sequences);
  }
  
  
  if ((length(sequences) == 0) || (nchar(prev_mer) > max.mer)) {
    return(file.cons);
  }
  for (base in bases) {
    new_mer = paste0(prev_mer, base);
    new_sequences = sequences[grep(new_mer, sequences)]
    if (length(new_sequences) > 0) {
      mer_size = nchar(new_mer);
      if (mer_size >= min.mer && mer_size <= max.mer) {
        label = paste0("mer",mer_size);
        if (!(label %in% names(file.cons))) {
          #cat("Creating ",label,"\n");
          file.cons[[label]] = file(label, "w");
        }
        #cat(new_mer,"\n");
        writeLines(new_mer, file.cons[[label]]);
        #flush(file.cons[[label]]);
      }
      
      #print(names(file.cons));
      file.cons = zcurve$writeMers(
        sequences = new_sequences,
        prev_mer = new_mer,
        bases = bases,
        min.mer = min.mer,
        max.mer = max.mer,
        file.cons = file.cons
      );
    }
  }
  #cat("Outside loop\n");
  #print(names(file.cons));
  if (prev_mer == "") { # First call, flush and close.
    for (file.con in names(file.cons)) {
      flush(file.cons[[file.con]]);
      close(file.cons[[file.con]]);
    }
  }
  return(file.cons);
}



