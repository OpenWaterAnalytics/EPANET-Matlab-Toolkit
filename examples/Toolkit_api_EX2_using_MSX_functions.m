
%Clear
clear; close('all'); clc;
start_toolkit;
% Create EPANET object using the INP file
inpname='example.inp'; %net2-cl2 example
%% apiMSX
d=epanet(inpname);
msxname = [inpname(1:end-4),'.msx'];
d.apiMSXMatlabSetup(d, msxname);


d.apiMSXsolveH(d.MSXLibEPANET);

d.apiMSXinit(0, d.MSXLibEPANET);

ss=1:d.LinkCount;%index link
uu=1:d.MSXSpeciesCount;
value.Quality = cell(1, length(ss));
% Obtain a hydraulic solution
d.apiMSXsolveH(d.MSXLibEPANET);
% Run a step-wise water quality analysis without saving
% RESULTS to file
d.apiMSXinit(0, d.MSXLibEPANET);
% Retrieve species concentration at node
k = 1; tleft = 1;
time_step = 360;
value.Time(k, :)=0;
            while(tleft>0)
                [t, tleft]= d.apiMSXstep(d.MSXLibEPANET);
                if t<time_step || t==time_step
                    i=1;
                     for nl=ss
                        g=1;
                        for j=uu
                            [obj.Errcode, value{i}(j)] = obj.apiMSXgetinitqual(1, i, j, obj.MSXLibEPANET);
                    if obj.Errcode, error(obj.getMSXError(obj.Errcode)); end
                            value.Quality{i}(k, g)=obj.getMSXLinkInitqualValue{(nl)}(j);
                            g=g+1;
                        end
                        i=i+1;
                    end
                else
                    i=1;
                    for nl=ss
                        g=1;
                        for j=uu
                            value.Quality{i}(k, g)=obj.getMSXSpeciesConcentration(1, (nl), j);%link code1
                            g=g+1;
                        end
                        i=i+1;
                    end
                end
                end
                if k>1
                    value.Time(k, :)=t;
                end
                k=k+1;
