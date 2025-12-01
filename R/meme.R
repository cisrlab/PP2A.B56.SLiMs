

have_ext_fimo<-function() {
  return(file.exists(itcpredictr:::meme$getBin("fimo")))
}


#https://www.r-bloggers.com/identifying-the-os-from-r/
get_os<-function() {
  sysinf <- Sys.info();
  os <- "Unknown"
  if (!is.null(sysinf)) {
    os <- sysinf['sysname']
    if (os == "Darwin") {
      os <- "osx";
    }
  } else {
    os <- .Platform$OS.type;
    if (grepl("^darwin", R.version$os)) {
      os <- "osx";
    }
    if (grepl("linux-gnu", R.version$os)) {
      os <- "linux";
    }
  }
  return(tolower(os))
}

meme=list();

#Finds the path of the binary,
#For windows, assumes that meme has been installed under ~/meme of the
#cygwin directory.
meme$getBin <- function(cmd_name) {
  if (get_os() == "windows") {
    username = Sys.info()[["user"]]
    cmd = paste0("C:/cygwin64/bin/bash.exe --login -c \"/home/",
      username,
      "/bin/", cmd_name, ".exe");
    return(cmd);
  } else {
    ret = system(paste0("which ",cmd_name), intern=TRUE)
    if (length(ret) == 0) {
      stop("Cannot find binary ",cmd_name)
    }
    return(ret)
  }
}

meme$getMeme <- function() {
  return(meme$getBin("meme"));
}

meme$getFimo <- function() {
  return(meme$getBin("fimo"));
}

meme$meme<-function(
  sequences,
  weights = c(),
  dir=paste0(getwd(),"/temp.meme"),
  mod = "zoops",
  minsites = NA,
  minw = NA,
  maxw = NA,
  nmotifs = 1,
  bg_file = "",
  debug=FALSE
  ) {
    if (debug) {cat("meme start\n");}

    if (debug) {print(class(sequences));}
    if (class(sequences) == "character") {
      library(Biostrings);
      sequences = BStringSet(sequences);
      names(sequences) = paste0("seq_",1:length(sequences));
    }

    fasta_file = tempfile();
    if (debug) {
      fasta_file = paste0(getwd(),"/temp.fasta");
    }
    unlink(fasta_file);
    writeXStringSet(sequences, fasta_file);

    if (length(weights)==length(sequences)) {
      weights = weights/max(weights);
      cat("usings weights");
      wline = ">WEIGHTS";
      for (idx in 1:length(weights)) {
        wline = paste(wline, weights[idx]);
      }
      con = file(fasta_file, "r");
      lines = readLines(con);
      close(con);
      lines = c(wline, lines);
      writeLines(lines, fasta_file);
    }

    meme_cmd = meme$getMeme();

    if (get_os() == "windows") {
      meme_cmd = paste(meme_cmd, "-oc", cygwinpath(dir), sep = " ");
    } else {
      meme_cmd = paste(meme_cmd, "-oc", paste0("\"", dir, "\""), sep = " ");
    }

    meme_cmd = paste(meme_cmd, "-nmotifs", nmotifs, "-protein", sep = " ");
    meme_cmd = paste(meme_cmd, "-mod", mod, sep = " ");
    if (!is.na(minsites)) {
      meme_cmd = paste(meme_cmd, "-minsites", minsites, sep = " ");
    }
    if (!is.na(minw)) {
      meme_cmd = paste(meme_cmd, "-minw", minw);
    }
    if (!is.na(maxw)) {
      meme_cmd = paste(meme_cmd, "-maxw", maxw);
    }
    if (bg_file != "") {
      meme_cmd = paste(meme_cmd, "-bfile", bg_file, sep = " ");
    }

    if (debug) {
      meme_cmd = paste(meme_cmd, "-V", sep= " ");
    }

    if (get_os() == "windows") {
      meme_cmd = paste(meme_cmd, cygwinpath(fasta_file), sep=" ");
    } else {
      meme_cmd = paste(meme_cmd, fasta_file, sep=" ");
    }
    if (debug) {
      cat(meme_cmd, "\n");
    }

    if (get_os() == "windows") {
      meme_cmd = paste0(meme_cmd,"\"");
    }


    unlink(dir, recursive=TRUE);
    ret = system(meme_cmd);

    file.copy(fasta_file, paste0(dir, "/seqs.fasta"));
    if (!debug) {
      unlink(fasta_file);
    }
    ans = paste0(dir, "/", "meme.xml");
    attr(ans, "png") = paste0(dir, "/", "logo1.png");
    attr(ans, "eps") = paste0(dir, "/", "logo1.eps");


    if (get_os() == "windows") {
      convert_cmd = paste0("C:/cygwin64/bin/bash --login -c \"/usr/bin/convert", " ",
        sub("C:/","/cygdrive/c/", attr(ans, "eps")), " ",
        sub("C:/","/cygdrive/c/", attr(ans, "png")), "\"");
      cat("cmd:",convert_cmd,"\n");
      ret = system(convert_cmd);

    }


    return(ans);

}

meme$fasta_get_markov<-function(sequence_file, background_file, order = 0) {

  meme_cmd = "fasta-get-markov";
  meme_cmd = paste(meme_cmd, "-m", order, sep=" ", collapse=" ");
  meme_cmd = paste(meme_cmd, sequence_file, background_file, sep=" ", collapse=" ");

  cat("cmd:", meme_cmd, "\n");
  ret = system(meme_cmd);

}


#--max-stored-scores <int> (default=100000)
meme$fimo<-function(
  seqs,
  meme_xml,
  bg_file="",
  dir=paste0(getwd(),"/fimo_output"),
  no_qvalue = FALSE,
  thresh=0.5,
  max_stored_scores = 100000,
  debug = FALSE,
  use_internal = !have_ext_fimo()
  ) {

  if (debug) {cat("meme$fimo: start\n");}
  library(Biostrings);
  if (class(seqs) == "character") {
    sequences = BStringSet(seqs);
    names(sequences) = paste0("seq_",1:length(sequences));
  } else {
    sequences = seqs;
  }
  if (debug) {
    file.path = "fimo.fasta"
  } else {
    file.path = tempfile();
  }
  unlink(file.path);
  writeXStringSet(sequences, file.path);

  fxn <- meme$fimoFasta
  if (use_internal) {
      message("Using internal meme")
      fxn <- itcpredictr:::fimoF
  }


  fimo.results = tryCatch(
    expr = fxn(
      file.path,
      meme_xml = meme_xml,
      bg_file = bg_file,
      dir = dir,
      thresh = thresh,
      no_qvalue = no_qvalue,
      max_stored_scores = max_stored_scores,
      debug = debug
    ),
    finally = function() {
        message("Finally, ", debug)
        if (!debug){unlink(file.path, recursive = TRUE)}
    }
  );
  return(fimo.results);
}

cygwinpath<-function(window_path, quotes=TRUE) {
  if (get_os() == "windows") {
    ans = sub("C:/", "/cygdrive/c/", window_path);
    ans = sub("C:\\\\", "/cygdrive/c/", ans);
    ans = sub("E:/", "/cygdrive/e/", ans);
    ans = sub("E:\\\\", "/cygdrive/e/", ans);

    ans = gsub("\\\\","/", ans);
    if (quotes) {
      ans = paste0("'",ans,"'");
    }
    return(ans);
  } else {
    return(window_path);
  }
}

#
meme$fimoFasta<-function(
  fasta,
  meme_xml,
  bg_file="",
  dir = "fimo_output",
  no_qvalue = FALSE,
  thresh=0.5,
  max_stored_scores = 100000,
  debug = FALSE
  ) {

  if (debug) {cat("meme$fimoFasta: start\n");}

  fimo_cmd = meme$getFimo();

  fimo_args=paste("--norc", "--thresh",thresh, sep = " ");
  if (get_os() == "windows") {
    fimo_args = paste(fimo_args, "--oc", cygwinpath(dir), sep = " ");
  } else {
    fimo_args = paste(fimo_args, "--oc", paste0("\"", dir, "\""), sep = " ");
  }
  if (bg_file != "") {
     fimo_args = paste(fimo_args, "--bgfile", bg_file, sep = " ");
  }
  if (no_qvalue) {
    fimo_args = paste(fimo_args, "--no-qvalue");
  }

  fimo_args = paste(fimo_args, "--max-stored-scores", max_stored_scores)
  #print(fimo_args);
  if (get_os() == "windows") {
    fimo_args = paste(fimo_args, cygwinpath(meme_xml), cygwinpath(fasta), sep= " ");
    fimo_args = paste0(fimo_args, "\"");
  } else {
    fimo_args = paste(fimo_args, meme_xml, fasta,sep=" ");
  }
  if (debug) {cat(fimo_cmd," ", fimo_args, "\n");}
  unlink(dir);

  if (dir.exists(dir)) {
    unlink(dir, recursive=TRUE);
  }
  dir.create(dir);

  ret = system(paste0(fimo_cmd, " ", fimo_args), ignore.stderr = !debug, ignore.stdout = !debug)
  if (ret != 0) {
    stop("Error running fimo:", fimo_cmd, " ",fimo_args, " ", ret);
  }
  fimo_txt = data.frame();
  if (file.exists(paste0(dir,"/fimo.txt"))) {
    fimo_txt = readFimoTxt(paste0(dir,"/fimo.txt"), debug=debug);
  } else if (file.exists(paste0(dir,"/fimo.tsv"))) {
    fimo_txt = readFimoTxt(paste0(dir,"/fimo.tsv"), debug=debug);
  } else {
    stop("Cannot find fimo results");
  }
  if(debug){print(head(fimo_txt))}
  if (no_qvalue) {
    fimo_txt = fimo_txt[,colnames(fimo_txt) != "q.value"];
  }

  if (dir.exists(dir) && !debug) {
      unlink(dir, recursive=TRUE);
  }

  return(fimo_txt);
}

#
# Attempts to read in the fimo.txt file
#
readFimoTxt<-function(
  fimo_txt,
  debug=FALSE
  ) {

  if (debug) { cat("readFimoTxt: start\n"); }
  if (debug) { cat("readFimoTxt: reading ",fimo_txt,"\n"); }
  con = file(fimo_txt,"r");
  lines = readLines(con);
  close(con);
  ans <- readFimoLines(lines, debug=debug)
  return(ans)
}

readFimoLines<-function(lines, debug=FALSE) {
  if (debug) { print("raw file"); print(head(lines)); }

  lines = lines[lines != ""];
  lines = lines[grep("# FIMO",lines,invert=TRUE)];
  lines = lines[grep("# fimo",lines,invert=TRUE)];
  lines = lines[grep("# The format",lines,invert=TRUE)];
  lines = lines[grep("# /",lines,invert=TRUE)];
  if (debug) { print("filtered lines"); print(head(lines));}

  tokens = strsplit(lines,"\t");
#  if (debug) { print("tokens"); print(tokens);}
  if (length(tokens) == 0) {return (data.frame());}
  ncol = length(tokens[[1]]);

  ans.df = t(matrix(unlist(tokens),nrow=ncol));

  has_header <- ans.df[1,1] == "#pattern name"

  # First time, fimo internal call, header exists.  Next call,
  # Header is not printed.  This is due to an internal static
  # variable within fimo-output.c
  if (has_header) {
    header = ans.df[1,];
    if (debug) { print("header"); print(header); }
    if (debug) { print("ans.df"); print(head(ans.df)); }
    ans.df = ans.df[-1,];
  } else {
    header_cols <- c(
      "#pattern name",
      "sequence name",
      "start",
      "stop",
      "strand",
      "score",
      "p-value",
      "q-value",
      "matched sequence"
    )
    header <- header_cols
  }

  if (debug) { print("ans.df2"); print(head(ans.df)); }
  #print(ans.df);
  if (debug) {
    message("exists:",exists("ans.df"))
    message("inherits(ans.df,\"matrix\"):", inherits(ans.df,"matrix"))
    message("class(ans.df):",class(ans.df))
    message("is.null:",is.null(ans.df))
    message("length(ans.df):",length(ans.df))
    message("nrow(ans.df):",nrow(ans.df))
    message("is.na(nrow(ans.df)):", is.na(nrow(ans.df)))
  }
  if (!exists("ans.df") ||
      is.null(ans.df) ||
      length(ans.df) == 0 ||
      !inherits(ans.df,"matrix") ||
      nrow(ans.df) == 0) {
    return(data.frame());
  }
  if (nrow(ans.df) == 1) {
    ans.df = as.data.frame(t(ans.df), stringsAsFactors=FALSE);
  } else {
    ans.df = as.data.frame(ans.df, stringsAsFactors=FALSE);
  }
  if (debug) { print("ans.df3"); print(head(ans.df)); }
  colnames(ans.df) = header;

  if (debug) { print("ans.df4"); print(head(ans.df)); }

  colnames(ans.df)[colnames(ans.df) == "# motif_id"] = "motif.id";
  colnames(ans.df)[colnames(ans.df) == "motif_id"] = "motif.id";
  colnames(ans.df)[colnames(ans.df) == "#pattern name"] = "motif.id";
  colnames(ans.df)[colnames(ans.df) == "motif_alt_id"] = "motif.alt.id";
  colnames(ans.df)[colnames(ans.df) == "matched_sequence"] = "motif.sequence";
  colnames(ans.df)[colnames(ans.df) == "matched sequence"] = "motif.sequence";
  colnames(ans.df)[colnames(ans.df) == "matched.sequence"] = "motif.sequence";
  colnames(ans.df)[colnames(ans.df) == "p-value"] = "p.value";
  colnames(ans.df)[colnames(ans.df) == "q-value"] = "q.value";
  colnames(ans.df)[colnames(ans.df) == "start"] = "motif.start";
  colnames(ans.df)[colnames(ans.df) == "stop"]  = "motif.end";
  colnames(ans.df)[colnames(ans.df) == "sequence_name"] = "sequence.name";
  colnames(ans.df)[colnames(ans.df) == "sequence name"] = "sequence.name";

  ans.df$p.value = as.numeric(ans.df$p.value);
  ans.df$q.value = as.numeric(ans.df$q.value);
  ans.df$score = as.numeric(ans.df$score);
  ans.df$motif.start = as.integer(ans.df$motif.start);
  ans.df$motif.end = as.integer(ans.df$motif.end);
  ans.df = ans.df[!is.na(ans.df$motif.start),];
  if (debug) { cat("readFimoTxt: stop\n");}
  if (debug) {
    attr(ans.df, "lines") = lines;
    #attr(ans.df, "tokens") = tokens;
  }

  ans.df <- ans.df[order(ans.df$p.value, decreasing=FALSE),]


  return(ans.df);


}

meme$fimoSeqChr<-function(seq_char, meme_xml, bg_file="", no_qvalue = FALSE, thresh = 0.5, debug = FALSE) {

  library(Biostrings)

  x = BStringSet(seq_char)
  names(x)[1] = "temp";
  writeXStringSet(x, "temp.fasta");
  a = meme$fimoFasta("temp.fasta", meme_xml, bg_file, no_qvalue = no_qvalue, thresh=thresh, debug = debug);
  #a$sequence = seq_char;
  return(a);

}

meme$iupac2meme<-function(iupac_seqs, meme_out, numseqs=20, debug=FALSE) {
  cmd = paste(meme$getBin("iupac2meme"), "-protein");
  cmd = paste(cmd, "-numseqs", numseqs);
  cmd = paste(cmd, paste(iupac_seqs, sep=" ", collapse=" "));
  cmd = paste(cmd, "> ",file.path);

  if (debug) {
    cat("cmd:",cmd,"\n");
  }

  ret = system(cmd);

  if (debug) {
    cat("ret:", ret, "\n");
  }

  ans = file.path;
  attr(ans, "ret") = ret;
  attr(ans, "cmd") = cmd;
  return(ans);

}

meme$meme2png<-function(meme_xml, motif_id, output_path="logo1.png", debug=FALSE) {
  # TODO fix this to work on both windows, mac, and linux with temporary directories.

  cmd = paste(meme$getBin("meme2images"));

  tdir = "meme_png";
  if (!missing(motif_id)) {
    cmd = paste(cmd, "-motif", motif_id);
  }
  cmd = paste(cmd, "-png");
  cmd = paste(cmd, cygwinpath(meme_xml));
  cmd = paste0(cmd, " '", tdir, "'");

  if (debug) {
    cat("cmd:",cmd,"\n");
  }
  ret = system(cmd);

  if (ret != 0) {
    stop("Error running meme2images! ", ret);
  }

  print( list.files(paste0(tdir,"/")));

  image_list = list.files(paste0(tdir, "/"));

  src_path = cygwinpath(paste0(tdir, "/" , image_list), quotes=FALSE);

  print(src_path);
  for (file_idx in 1:1) {
    ret = file.copy(src_path[file_idx],output_path)
    if (!ret) {
      warning("Error copying ", src_path, " to ", output_path , " ", ret);
    }
  }
  unlink(tdir, recursive=TRUE);

}

meme$meme2meme<-function(meme_xmls, output_path="meme.xml", numbers=FALSE, debug=FALSE) {

  if (length(meme_xmls) < 1) {
    stop("Need at least one meme file!");
  }
  cmd = meme$getBin("meme2meme");

  args = "";

  if (numbers) {
    args = paste(args, "-numbers");
  }

  for (idx in 1:length(meme_xmls)) {
    args = paste(args, cygwinpath(meme_xmls[idx]));
  }

  temp_output = tempfile();


  if (debug) {
    cat("cmd:", cmd, "\n");
    cat("args:", args, "\n");
    cat("output_path:", cygwinpath(temp_output), "\n");
  }


#  ret = system2(paste(cmd,args,"\""), "", stdout = cygwinpath(temp_output));
  cat("cmd:",paste(cmd,args,"\""),"\n");
  if (get_os() == "windows") {
    ret = system(paste(cmd,args,">", cygwinpath(temp_output),"\""));
  } else {
    ret = system(paste(cmd,args,">", temp_output))
  }
  if (ret != 0) {
    stop("Error running meme2meme! ",ret);
  }
  file.copy(temp_output, output_path);
  return(output_path);
}


meme$testMeme<-function(debug=FALSE) {
  seqs = c(
  "KPPISPPKTPVPQASS",
  "KPPIAPKPVIPQLPTT",
  "PAKPPVSPKPVLTQPV",
  "YHSPVSRLPPSGSKPV",
  "QGRPPKVPVREVCPVT",
  "TVPDSCFPATKPPLPH",
  "EPIIPFRETITKPPKV",
  "MAADVSVTHRPPLSPK",
  "TTFEVPVSVHTRPPMK",
  "KPPIAAHASRSAESKT",
  "VTETVTTRLTSLPPKG",
  "AIISSSEKLLAKKPPS",
  "KPPMNGYQKFSQELLS",
  "PLPEPGYFTKPPIAAH",
  "TKPPLPHAACHSCSED",
  "GSKPPLPTSGPAQVPT",
  "VDTKPPVAHTNHILKL",
  "LPPPSPAKSDTLIVPS",
  "VPSEKPPMMPQAQPAI",
  "YDRKPPSGFKPLTPPT",
  "SSQDPANLFSLPPLPP",
  "LPPKPARPPAVSPALT",
  "RPPLPESLSTFPNTIS",
  "KPPLPMGSQVLQIRPN",
  "PQSRPPIPRTQPQPEP",
  "LPPSHGSSSGHPSKPY",
  "EPVKPAVLAKPPVCPA",
  "SPAVSQKPPFQSGMHY",
  "SKKPPPVPALPSKLPT",
  "LPPTTSSSKKPIPTLA",
  "IPPKVPQRTTSISPAL",
  "EPLSSKPPLPRKPLLQ",
  "DGHTFLLEKPPVPPKP"
  );

  meme.results = meme$meme(seqs, debug=debug);
  fimo.results = meme_fimo(seqs, meme.results, debug=debug);
  fimo.results.internal = meme_fimo(seqs, meme.results, debug=debug, use_internal = TRUE)
  ans = list();
  ans$meme.results = meme.results
  ans$fimo.results = fimo.results
  ans$fimo.results2 = fimo.results.internal
  return(ans);
}

get_motif_ids<-function(meme_model) {

  motif_lines = meme_get_lines(meme_model);
  motif_lines = motif_lines[grep("MOTIF",motif_lines)];

  motif_df = data.frame(
    motif.idx = 1:length(motif_lines),
    motif.id  = rep(NA, length(motif_lines)),
    motif.alt.id = rep(NA, length(motif_lines)),
    stringsAsFactors=FALSE
  );

  for(idx in 1:length(motif_lines)) {
    tokens = strsplit(motif_lines[idx]," ")[[1]];
    print(tokens);
    if (length(tokens) > 1) {
      motif_df$motif.id[idx] = tokens[2];
    }
    if (length(tokens) > 2) {
      motif_df$motif.alt.id[idx] = tokens[3];
    }
  }

  return(motif_df);
}

meme_nmotifs<-function(meme_model) {

  return(nrow(get_motif_ids(meme_model)));

}


meme_get_lines<-function(meme_model) {
  cmd = meme$getBin("meme-get-motif");

  args = paste("-all", cygwinpath(meme_model));
  print(cmd);
  print(args);
  temp_file = tempfile();
  cmd = paste(cmd,args, ">", cygwinpath(temp_file));
  if (get_os() == "windows") {
    cmd = paste(cmd, "\"");
  }
  print(cmd);

  ret = system(cmd);

  if (ret == 0) {
    file.con = file(temp_file,"r");
    lines = readLines(file.con);
    close(file.con);
    return(lines);
  } else {
    stop("Error :", ret);
  }
}

meme_motif_change_ids<-function(meme_model_in, meme_model_out, motif_indices, motif_ids, motif_alt_ids) {

  motif_lines = meme_get_lines(meme_model_in);
  motif_line_indices = grep("MOTIF", motif_lines);
  for (motif_idx_idx in 1:length(motif_indices)) {
    motif_idx = motif_indices[motif_idx_idx];
    motif_line_idx = motif_line_indices[motif_idx];
    motif_lines[motif_line_idx] = paste("MOTIF", motif_ids[motif_idx_idx], motif_alt_ids[motif_idx_idx], sep=" ", collapse=" ");
  }
  file.con = file(meme_model_out, "w");
  writeLines(con = file.con, text = motif_lines)
  close(file.con);

}



