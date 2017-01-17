void printSyms(std::vector< std::vector<size_t> >& syms) {
	int ngid = syms.size();

	for(int gid=0; gid<ngid; ++gid) {
		for(int i=0; i<syms[gid].size(); ++i) {
			cout << syms[gid][i] << " " ;
		}
		cout << endl;
	}
}