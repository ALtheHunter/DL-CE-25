addpath('../main')
addpath('../stateEvo')
addpath('../classification')

% handle random seed
if verLessThan('matlab','7.14')
  defaultStream = RandStream.getDefaultStream;
else
  defaultStream = RandStream.getGlobalStream;
end
if 1 % new RANDOM trial
  savedState = defaultStream.State;
  save random_state.mat savedState;
else % repeat last trial
  load random_state.mat
end
defaultStream.State = savedState;

% simulation parameters
L = 100;  % # of measurement vectors [100]
switch 1
 case 1
  likeType = 'Probit'; % in {'AWGN','Probit'}
  N = 1024; % signal dimension [1024]
  del = 4.0; % measurement rate M/N [4.0]
  beta = 1/32; % sparsity rate K/N [1/32]
%   N = 1024; % signal dimension [1024]
%   del = 0.7; % measurement rate M/N [4.0]
%   beta =0.01; % sparsity rate K/N [1/32]
 case 2
  likeType = 'AWGN'; % in {'AWGN','Probit'}
  N = 512; % signal dimension [512]
  del = 0.5; % measurement rate M/N [0.5]
  beta = 1/32; % sparsity rate K/N [0.1]
end
SNRdB = 40; % [40]
svType = 'cond_num'; % in {'cond_num','spread','low_rank'}
cond_num = 1; % condition number
spread = 1; % amount to spread singular values (=1 means iid Gaussian A, =0 means frame)
low_rank = round(min(N,round(del*N))/8);
UType = 'Haar'; % in {'DFT','DCT','DHT','DHTrice','Haar','I'}
VType = 'Haar'; % in {'DFT','DCT','DHT','DHTrice','Haar','I'}
shuffle = true; % shuffle rows of V' ?
randsign = true; % randomly sign-flip columns of V' ?
isCmplx = false; % simulate complex-valued case?
plot_traj = true; % plot trajectory of each column?
plot_sig = true; % plot trajectory of each column?
runOracle = true;  % calculate support oracle?
runGAMP = true;  % run GAMP?

% algorithmic parameters
maxIt = 50; % max iterations for VAMP
tol = min(1e-3,max(1e-6,10^(-SNRdB/10))); % stopping tolerance for VAMP
damp = 0.9; % damping parameter
denoiser = 'BG'; % in {'BG','DMM','MAPLaplace'}
learnPrior = false; % automatically tune the denoiser?
learnLike = false; % automatically tune the likelihood?

% other defaults
fixed_K = true; % used fixed sparsity K=E{K}=round(rho*M)?
Afro2 = N; % squared Frobenius norm of matrix
xvar0 = 1; % prior variance of x elements
xmean1 = 0; % prior mean of non-zero x coefs

% setup
M = round(del*N);
xvar1 = xvar0/beta; % prior variance of non-zero x coefs
wvar = (Afro2/M)*10^(-SNRdB/10)*beta*(abs(xmean1)^2 + xvar1); 
if (strcmp(UType,'DFT')||strcmp(VType,'DFT'))&&(~isCmplx)
  warning('setting isCmplx=true since complex-valued matrix')
  isCmplx = true;
elseif isCmplx&&(~strcmp(UType,'DFT'))&&(~strcmp(VType,'DFT'))
  warning('setting isCmplx=false since real-valued matrix')
  isCmplx = false;
end

% generate signal 
x = zeros(N,L); 
for l=1:L
  if fixed_K
    supp = randperm(N,round(beta*N));
  else
    supp = find(rand(N,1)<beta);
  end
  K = length(supp);
  %supp = 1:K; display('block support'); % for testing
  if isCmplx
    x(supp,l) = xmean1 + sqrt(0.5*xvar1)*randn(K,2)*[1;1j];% complex Gaussian
  else
    x(supp,l) = xmean1 + sqrt(xvar1)*randn(K,1);% real Gaussian
  end
end
%x = abs(x); display('positive x')
%x =x(:,1)*ones(1,L); display('repeated x')

% generate noise
if isCmplx
  w = sqrt(0.5*wvar)*(randn(M,L) + 1j*randn(M,L));
else
  w = sqrt(wvar)*randn(M,L);
end

% generate linear transform
switch svType
  case 'spread', svParam = spread;
  case 'cond_num', svParam = cond_num;
  case 'low_rank', svParam = low_rank;
end
mat = genMatSVD(M,N,UType,svType,svParam,VType,...
                'isCmplx',isCmplx,'Afro2',Afro2,...
                'shuffle',shuffle,'randsign',randsign,...
                'fxnHandles',true);
A = mat.fxnA; Ah = mat.fxnAh;
U = mat.fxnU; Uh = mat.fxnUh;
V = mat.fxnV; Vh = mat.fxnVh;
d = mat.s.^2; 

% crease noisy observations
%z = A*x; 
z = A(x); 
SNRdB_test = 20*log10(norm(z(:))/norm(w(:)));
switch likeType
  case 'AWGN'
    y = z + w;
  case 'Probit'
    if isCmplx
      error('Set isCmplx=false for Probit likelihood')
    else
      y = ((z+w)>0);
    end
end

% support-oracle performance bound for AWGN case
if strcmp(likeType,'AWGN')&&runOracle
  I = speye(N);
  x0 = zeros(N,L);
  oracleNMSEdB = nan(1,L);
  for l=1:L
    supp = find(x(:,l)~=0);
    try
      A0 = A(I(:,supp)); % fast but not compatible with all A(.)
    catch
      K = length(supp); % slow but compatible with all A(.)
      A0 = zeros(M,K);
      for k=1:K, A0(:,k) = A([zeros(supp(k)-1,1);1;zeros(N-supp(k),1)]); end
    end
    a0 = A0*ones(length(supp),1);
    x0(supp,l) = xmean1 + (A0'*A0+(wvar/xvar1)*eye(length(supp)))\(A0'*(y(:,l)-a0*xmean1));
    oracleNMSEdB(l) = 20*log10(norm(x0(:,l)-x(:,l))/norm(x(:,l)));
  end
end

% establish input denoiser
switch denoiser
case 'BG'
  if learnPrior
    betaInit = 1/N; 
    xvar0init = xvar0;
    xvar1init = xvar0init/betaInit;
    tuneDim = 'joint';
    if isCmplx
      EstimIn = SparseScaEstim(CAwgnEstimIn(0,xvar1init,0,'autoTune',true,'mean0Tune',false,'tuneDim',tuneDim),betaInit,0,'autoTune',true,'tuneDim',tuneDim);
    else
      EstimIn = SparseScaEstim(AwgnEstimIn(0,xvar1init,0,'autoTune',true,'mean0Tune',false,'tuneDim',tuneDim),betaInit,0,'autoTune',true,'tuneDim',tuneDim);
    end
  else
    if isCmplx
      EstimIn = SparseScaEstim(CAwgnEstimIn(xmean1,xvar1),beta);
    else
      EstimIn = SparseScaEstim(AwgnEstimIn(xmean1,xvar1),beta);
    end
  end
case 'DMM'
  alpha = 1.5;
  debias = false; 
  EstimIn = SoftThreshDMMEstimIn(alpha,'debias',debias);
  if learnPrior, 
    warning('learnPrior not implemented for SoftThreshDMM'); 
  end;
case 'MAPLaplace'
  lam = 1/sqrt(wvar);
  if learnPrior,
    EstimIn = SoftThreshEstimIn(lam,0,'autoTune',true,'counter',10) 
  else
    EstimIn = SoftThreshEstimIn(lam);
  end
otherwise
  error('unknown denoiser')
end

% establish likelihood
tuneDim = 'joint';
switch likeType
  case 'AWGN'
    if learnLike
      wvarInit = 0.01 * norm(y,'fro')^2 / (M*L); % SNR ~= -20 dB
      if isCmplx,
        EstimOut = CAwgnEstimOut(y,wvarInit,false,'autoTune',true,...
                        'tuneMethod','EM','tuneDamp',1,'tuneDim',tuneDim);
      else
        EstimOut = AwgnEstimOut(y,wvarInit,false,'autoTune',true,...
                        'tuneMethod','EM','tuneDamp',1,'tuneDim',tuneDim);
      end
    else
      if isCmplx,
        EstimOut = CAwgnEstimOut(y,wvar);
      else
        EstimOut = AwgnEstimOut(y,wvar);
      end
    end
  case 'Probit'
    if learnLike, 
      wvarInit = 1e-10; % better choice?
      EstimOut = ProbitEstimOut(y,0,wvarInit,false,'autoTune',true,...
                        'tuneMethod','EM','tuneDim',tuneDim);
    else
      EstimOut = ProbitEstimOut(y,0,wvar);
    end
end

% setup VAMP
vampOpt = VampGlmOpt;
vampOpt.nitMax = maxIt;
vampOpt.tol = tol;
vampOpt.damp = damp;
vampOpt.verbose = false;
vampOpt.fxnErr1 = @(x1,z1) 10*log10( sum(abs(x1-x).^2,1)./sum(abs(x).^2,1) );
vampOpt.fxnErr2 = @(x1,z1) 10*log10( sum(abs(...
        bsxfun(@times, x1, sum(conj(x1).*x,1)./sum(abs(x1).^2,1)) - x...
        ).^2,1)./sum(abs(x).^2,1) );
vampOpt.Ah = Ah; vampOpt.d = d; vampOpt.N = N;
vampOpt.U = U; vampOpt.Uh = Uh; vampOpt.V = V; vampOpt.Vh = Vh;
if strcmp(denoiser,'DMM') % can't initialize at r1=0 !
  vampOpt.gam1xinit = 1;
  vampOpt.r1init = randn(N,1);
end
if isCmplx, 
  vampOpt.r1init(1) = eps*1i; % so that SparseScaEstim knows it's complex-valued
end

% run VAMP
[x1,vampEstFin] = VampGlmEst(EstimIn,EstimOut,A,vampOpt);
vampNMSEdB_ = vampEstFin.err1;
vampNMSEdB_debiased_ = vampEstFin.err2;
vampNit = vampEstFin.nit;

% run VAMP state evolution
estInAvg = EstimInAvg(EstimIn,x);
phat = sqrt(xvar0*Afro2/M)*randn(M,L);
switch likeType
  case 'AWGN'
    clear estOutAvg
    estOutAvg.mse = @(pvar) deal( 1/(1/wvar+1/pvar), 1/(1/wvar+1/pvar) ); 
    %estOutAvg = EstimOutAvg2(EstimOut,z); % monte-carlo given z
    %estOutAvg = EstimOutAvg(likeType,wvar,phat); % monte-carlo given phat
  case 'Probit'
    %estOutAvg = EstimOutAvg2(EstimOut,z); % monte-carlo given z
    estOutAvg = EstimOutAvg(likeType,wvar,phat); % monte-carlo given phat
end
vampSeNMSE = VampGlmSE(estInAvg,estOutAvg,d,N,M/N,vampNit)/(beta*(abs(xmean1)^2+xvar1));

% print learned parameters
if learnPrior
  beta
  betaEstimate = EstimIn.p1(1,:)
  xvar1
  xvar1estimate = EstimIn.estim1.var0(1,:)
end
if learnLike
  wvar
  switch likeType
    case 'AWGN'
      wvarEstimate = EstimOut.wvar(1) 
    case 'Probit'
      wvarEstimate = EstimOut.Var(1) 
  end
end

% setup and run GAMP
gampNit = 0;
if runGAMP
  %Agamp = MatrixLinTrans(A);
  Agamp = FxnhandleLinTrans(M,N,A,Ah,Afro2/(M*N));
  optGAMP = GampOpt('legacyOut',false,'step',0.25,'stepMax',0.25,'adaptStep',false,'uniformVariance',true,'tol',tol,'nit',500);
  if learnPrior
    % reset these values
    EstimIn.p1 = betaInit;
    EstimIn.estim1.var0 = xvar1init;
    EstimIn.estim1.mean0 = 0;
  end
  if learnLike
    % reset these values
    switch likeType
      case 'AWGN'
        EstimOut.wvar = wvarInit;
        EstimOut.tuneMethod = 'ML';
        EstimOut.tuneDim = 'col'; % seems to be important
      case 'Probit'
        EstimOut.Var = wvarInit;
        EstimOut.tuneMethod = 'ML';
        warning('NEED TO SET tuneDim=col')
    end
  end
  tstart = tic;
  [gampEstFin,optGampFin,gampEstHist] = gampEst(EstimIn,EstimOut,Agamp,optGAMP);
  time_gamp = toc(tstart);
  gampNit = gampEstFin.nit;
  gampNMSEdB_ = nan(L,gampNit);
  gampNMSEdB_debiased_ = nan(L,gampNit);
  for l=1:L
    xhat_ = gampEstHist.xhat((l-1)*N+[1:N],:);
    gampNMSEdB_(l,:) = 10*log10(sum(abs(xhat_-x(:,l)*ones(1,gampNit)).^2,1)/norm(x(:,l))^2);
    gain_ = conj(x(:,l)'*xhat_)./sum(abs(xhat_).^2,1);
    gampNMSEdB_debiased_(l,:) = 10*log10(sum(abs( bsxfun(@times,xhat_,gain_)-x(:,l)*ones(1,gampNit)).^2,1)/norm(x(:,l))^2);
  end
  %figure(4); clf; gampShowHist(gampEstHist,optGampFin,x); % debug GAMP
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
if plot_traj
figure(3); clf;
subplot(211) 
  % plot VAMP
  handy = plot(1:vampNit-1,vampNMSEdB_(:,2:end).','b.-'); % 1st iteration is trivial
  set(handy(1),'Displayname','VAMP') 
  % plot GAMP
  if runGAMP
    hold on; 
      handy = [handy, plot(1:gampNit-1,gampNMSEdB_(:,2:end),'r.-')];
      set(handy(1,end),'Displayname','GAMP') 
    hold off; 
  end
  % plot support oracle
  if strcmp(likeType,'AWGN')&&runOracle
    ax = gca; ax.ColorOrderIndex = 1; % use same colors
    hold on; 
      handy = [handy, plot([1;max(vampNit,gampNit)],[1;1]*oracleNMSEdB,'k--')]; 
      set(handy(1,end),'Displayname','oracle') 
    hold off; 
  end
  legend(handy(1,:))
  ylabel('NMSE [dB]')
  xlabel('iterations')
  grid on

subplot(212) 
  % plot VAMP
  handy = plot(1:vampNit-1,vampNMSEdB_debiased_(:,2:end).','b.-'); % 1st iteration trivial
  set(handy(1),'Displayname','VAMP') 
  % plot GAMP
  if runGAMP
    hold on; 
      handy = [handy, plot(1:gampNit-1,gampNMSEdB_debiased_(:,2:end),'r.-')];
      set(handy(1,end),'Displayname','GAMP') 
    hold off; 
  end
  % plot support oracle
  if strcmp(likeType,'AWGN')&&runOracle
    ax = gca; ax.ColorOrderIndex = 1; % use same colors
    hold on; 
      handy = [handy, plot([1;max(vampNit,gampNit)],[1;1]*oracleNMSEdB,'k--')]; 
      set(handy(1,end),'Displayname','oracle') 
    hold off; 
  end
  legend(handy(1,:))
  ylabel('debiased NMSE [dB]')
  xlabel('iterations')
  grid on
end % plot_traj

if plot_sig
  figure(2); clf;
  l = 1;
  if runGAMP, subplot(211); end;
    stem(x1(:,l))
    gain = sum(conj(x1(:,l)).*x(:,l),1)./sum(abs(x1(:,l)).^2);
    hold on; 
      stem(bsxfun(@times,x1(:,l),gain),'x'); 
      stem(x(:,l),'--'); 
    hold off;
    legend('VAMP','VAMP debiased','true')
    if L>1, title(['column ',num2str(l),' of ',num2str(L)]); end;
    xlabel('coefficient index')
    grid on;
  if runGAMP,
  subplot(212)
    xg = gampEstFin.xhat;
    stem(xg(:,l))
    gain = sum(conj(xg(:,l)).*x(:,l),1)./sum(abs(xg(:,l)).^2);
    hold on; 
      stem(bsxfun(@times,xg(:,l),gain),'x'); 
      stem(x(:,l),'--'); 
    hold off;
    legend('GAMP','GAMP debiased','true')
    xlabel('coefficient index')
    grid on;
  end
end

figure(1); clf;
  vampNMSE_avg = mean(10.^(vampNMSEdB_debiased_/10),1);
  if runGAMP, gampNMSE_avg = mean(10.^(gampNMSEdB_debiased_/10),1); end;
  plot(1:vampNit,vampNMSE_avg,'.-','Displayname','VAMP');
  set(gca,'YScale','log','XScale','log')
  hold on;
    if runGAMP, semilogx(1:gampNit,gampNMSE_avg,'.-','Displayname','GAMP'); end;
    semilogx(2:vampNit,vampSeNMSE(2:end),'k--','Displayname','VAMP-SE');
  hold off;
  axe = axis;
  axis([1,axe(2),10^floor(log10(min([vampNMSE_avg,vampSeNMSE]))),1])
  legend(gca,'show')
  grid on
  xlabel('iteration')
  ylabel('avg debiased NMSE')

return

figure(4); clf;
  vampNMSE_avg = mean(10.^(vampNMSEdB_/10),1);
  if runGAMP, gampNMSE_avg = mean(10.^(gampNMSEdB_/10),1); end;
  plot(1:vampNit,vampNMSE_avg,'.-','Displayname','VAMP');
  set(gca,'YScale','log','XScale','log')
  hold on;
    if runGAMP, semilogx(1:gampNit,gampNMSE_avg,'.-','Displayname','GAMP'); end;
    semilogx(2:vampNit,vampSeNMSE(2:end),'k--','Displayname','VAMP-SE');
  hold off;
  axe = axis;
  axis([1,axe(2),10^floor(log10(min([vampNMSE_avg,vampSeNMSE]))),1])
  legend(gca,'show')
  grid on
  xlabel('iteration')
  ylabel('avg NMSE')
