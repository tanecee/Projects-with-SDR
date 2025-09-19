function p = parametersOFDM(Nsym)
    p.Nfft=256; 
    p.Nsym=Nsym; % Sabit 2 yerine, fonksiyonun girdisi olan Nsym'yi kullan
    p.actScs=p.Nfft/2; 
    p.dataScs=p.actScs/2;
    p.pilotScs=p.actScs/2;
    p.dataInd= (1:2:2*p.dataScs)+ p.Nfft/4;
    p.pilotInd = (2:2:2*p.pilotScs) + p.Nfft/4;
    
    % DC (1-tabanlı 129) çıkarıldı
    p.pilotInd = p.pilotInd(p.pilotInd ~= (p.Nfft/2+1));
    p.dataInd = p.dataInd(p.dataInd ~= (p.Nfft/2+1));
    p.pilotScs = length(p.pilotInd);
    p.dataScs = length(p.dataInd);

    p.pilot1=repmat([1 ; -1],p.pilotScs/2,1);
    p.pilot2=repmat([-1 ; 1],p.pilotScs/2,1);
    p.sync=zadoffChuSeq(8,255);
    p.cpLength=p.Nfft/4; 
    p.wformLength=(p.Nfft+p.cpLength)*p.Nsym;
    p.M=4; 
    p.sample_rate=20e6;
end
