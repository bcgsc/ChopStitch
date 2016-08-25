


bool getSeq(std::ifstream &uFile, std::string &line) {
    bool good=false;
    std::string hline;
    line.clear();
    if(opt::fastq) {
        good=getline(uFile, hline);
        good=getline(uFile, line);
        good=getline(uFile, hline);
        good=getline(uFile, hline);
    }
    else {
        do {
            good=getline(uFile, hline);
            if(hline[0]=='>'&&!line.empty()) break;// !line.empty() for the first rec
            if(hline[0]!='>')line+=hline;
        } while(good);
        if(!good&&!line.empty())
            good=true;
    }
    return good;
}
