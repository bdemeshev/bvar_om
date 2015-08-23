% Outputs
% (1) Error bands or marginal densities of b's -- a stacked vector of free parameters in A0
% (2) Error bands or marginal densities of g's -- a stacked vector of free parameters in A+
% (3) Error bands (imstp-by-nvar-by-nvar) or marginal densities (a stacked vector) of impulse responses.
% (4) Error bands (nfyr or nfqm-by-nvar) or marginal densities (a stacked vector) of unconditional out-of-sample forecasts

%*** Label for impulse responses
xlab = {'P'
        'I_d'
        'I_k'
        'y'};
ylab = xlab;


rnum=6;   % number of rows for subplot
cnum=2;   % number of columns for subplot

Indxdeng = 0;  % index for density graphs.  1: enable (and idensable error bands)
               %                            0: disenable (and enable error bands)
ncom = 10;   % The number of combined intervals.  A large ncom implies a large size of the bin.
%           % Must set ncom so that ninv can be divided

load outdistcnt      %yhatincntHD990   %outyhatincnt
       % ndraws1 ndraws2 timeminutes bRange


%(1)------------------------
%    Error bands or densities of b's
%(1)------------------------
if 0  %Aband
   %
   kIndx=[1 2 4]; %16 18 22 23 25]; %1:sum(n0)];  % index for selected variables.  Length(kIndx)<=size(bcnt,2).
   lval = [];   % length(kIndx)-element vector of the lowest values on the axis;
         % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
   hval = [];   % length(kIndx)-element vector of the highest values on the axis;
         % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
   [bpdf,bpo,bprob] = fn_histpdfcnt(bcnt,ndraws2,bRange(:,1),bhbin,ninv,Indxdeng,ncom,kIndx,lval,hval);
            % the densities (histograms) of variables
   [berrl,berrh,berrl1,berrh1] = fn_demarw(bpo,bprob);
            % .90 and .68 error bands
end


%(2)------------------------
%    Error bands or densities of g's
%(2)------------------------
if 0  %Apband
   %
   kIndx=[1 2 4];  %1:sum(np)];  % index for selected variables.  Length(kIndx)<=size(bcnt,2).
   lval = [];   % length(kIndx)-element vector of the lowest values on the axis;
         % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
   hval = [];   % length(kIndx)-element vector of the highest values on the axis;
         % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
   [gpdf,gpo,gprob] = fn_histpdfcnt(gcnt,ndraws2,gRange(:,1),ghbin,ninv,Indxdeng,ncom,kIndx,lval,hval);
            % the densities (histograms) of variables
   [gerrl,gerrh,gerrl1,gerrh1] = fn_demarw(gpo,gprob);
            % .90 and .68 error bands
end


%(3)------------------------
%    Error bands or densities of impulse responses
%(3)------------------------
if (IndxImf & Imfband)
   %
   kIndx=[1 100 imstp*nvar^2];  % index for selected variables.  Length(kIndx)<=size(bcnt,2).
   lval = [];   % length(kIndx)-element vector of the lowest values on the axis;
         % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
   hval = [];   % length(kIndx)-element vector of the highest values on the axis;
         % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
   [imfpdf,imfpo,imfprob] = fn_histpdfcnt(imfcnt,ndraws2,imfRange(:,1),imfhbin,ninv,Indxdeng,ncom,kIndx,lval,hval);
            % ninv+2-by-imstp*nvar^2 matrix of the densities (histograms) of variables
   [imferrl,imferrh,imferrl1,imferrh1] = fn_demarw(imfpo,imfprob);
            % .90 and .68 error bands
   %*** plot impulse responses
   imferrl = reshape(imferrl,imstp,nvar^2);
   imferrh = reshape(imferrh,imstp,nvar^2);
   imferrl1 = reshape(imferrl1,imstp,nvar^2);
   imferrh1 = reshape(imferrh1,imstp,nvar^2);
   figure
   fn_imcerrgraph(imfhat,imferrl1,imferrh1,nvar,imstp,xlab,ylab,1,[6 12 24 36])  % .68 bands.
   subtitle('.68 Error Bands')
   figure
   fn_imcerrgraph(imfhat,imferrl,imferrh,nvar,imstp,xlab,ylab,1,[6 12 24 36])  % .90 bands.
   subtitle('.90 Error Bands')
   %
   %============ For debugging ==========
   %[imferrl1(:,4) imfhat(:,4) imfmean(:,4) imferrh1(:,4)]
   [imferrl(1:2,9) imferrl1(1:2,9) imfhat(1:2,9) imfmean(1:2,9) imferrh1(1:2,9) imferrh(1:2,9)]  % R response to an MS shock.
end


%(4)------------------------
%    Error bands or densities of unconditional out-of-sample forecasts
%(4)------------------------
if 0   %(IndxFore & Foreband)
   %
   if 1    % annual grow rates
      kIndx=[1 6 nfyr*nvar];  % index for selected variables.  Length(kIndx)<=size(bcnt,2).
      lval = [];   % length(kIndx)-element vector of the lowest values on the axis;
            % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
      hval = [];   % length(kIndx)-element vector of the highest values on the axis;
            % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
      [yrgupdf,yrgupo,yrguprob] = fn_histpdfcnt(yrgucnt,ndraws2,yrguRange(:,1),yrguhbin,ninv,Indxdeng,ncom,kIndx,lval,hval);
               % ninv+2-by-imstp*nvar^2 matrix of the densities (histograms) of variables
      [yrguerrl,yrguerrh,yrguerrl1,yrguerrh1] = fn_demarw(yrgupo,yrguprob);
               % .90 and .68 error bands
      %*** plot impulse responses
      yrguerrl = reshape(yrguerrl,nfyr,nvar);
      yrguerrh = reshape(yrguerrh,nfyr,nvar);
      yrguerrl1 = reshape(yrguerrl1,nfyr,nvar);
      yrguerrh1 = reshape(yrguerrh1,nfyr,nvar);

      keyindx=[1:nvar];
      conlab=['unconditional'];
      figure
      fyrdates = [fdates(1,1):fdates(end,1)]';  % column vector
      fyrdates = [fyrdates zeros(size(fyrdates))];
      fn_forerrgraph(yafyrghate,[fyrdates yrguerrl1],[fyrdates yrguerrh1],...
                      yact2yrge,keyindx,rnum,cnum,q_m,ylab,forelabel,conlab)
   end
   %
   if 1  % period-to-period in the last year growth (annualized)
      kIndx=[1 6 nfqm*nvar];  % index for selected variables.  Length(kIndx)<=size(bcnt,2).
      lval = [];   % length(kIndx)-element vector of the lowest values on the axis;
            % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
      hval = [];   % length(kIndx)-element vector of the highest values on the axis;
            % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
      [qmyrgupdf,qmyrgupo,qmyrguprob] = fn_histpdfcnt(qmyrgucnt,ndraws2,qmyrguRange(:,1),qmyrguhbin,ninv,Indxdeng,ncom,kIndx,lval,hval);
               % ninv+2-by-imstp*nvar^2 matrix of the densities (histograms) of variables
      [qmyrguerrl,qmyrguerrh,qmyrguerrl1,qmyrguerrh1] = fn_demarw(qmyrgupo,qmyrguprob);
               % .90 and .68 error bands
      %*** plot impulse responses
      qmyrguerrl = reshape(qmyrguerrl,nfqm,nvar);
      qmyrguerrh = reshape(qmyrguerrh,nfqm,nvar);
      qmyrguerrl1 = reshape(qmyrguerrl1,nfqm,nvar);
      qmyrguerrh1 = reshape(qmyrguerrh1,nfqm,nvar);

      keyindx=[1:nvar];
      conlab=['unconditional'];
      figure
      fn_forerrgraph(yafqmyghate,[fdates qmyrguerrl1],[fdates qmyrguerrh1],...
                          yact2qmyge,keyindx,rnum,cnum,q_m,ylab,forelabel,conlab)
   end
end
