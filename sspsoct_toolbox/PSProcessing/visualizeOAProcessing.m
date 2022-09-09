function fh = visualizeOAProcessing(str,oaSheathStack,fhin,oaxOffset)
%fh = visualizeOAProcessing(str,fhin) visualizes the processing parameters
%of the structure str obtained by calling readFromUint8 with the optic axis
%channel of the respective pstif file.

if isfield(str,'tissueRet')
    Nrows = 4;
else
    Nrows = 2;
end

if nargin<2 || isempty(oaSheathStack)
    oaSheathStack = [];
end

if nargin<3 || isempty(fhin)
    fh = figure;
else
    fh = fhin;
end

if nargin<4 || isempty(oaxOffset)
    oaxOffset = [];
end

figure(fh);
clf
subplot(Nrows,2,1)
plot(str.wcorr')
xlabel('Spectral bins')
ylabel('[rad]')
set(gca,'xlim',[1,size(str.wcorr,2)])
title('Rotation vectors for symmetrization')
legend('Q','U','V')
    
subplot(Nrows,2,2)
hold on
plot(str.errInit*100)
plot(str.errFinal*100)
xlabel('Spectral bins')
ylabel('% of V component')
set(gca,'xlim',[1,size(str.wcorr,2)],'ylim',[0,40],'box','on')
title('Signal energy in V')
legend('Initial','Final')
    
subplot(Nrows,2,3)
if isfield(str,'alignErrInit')
    plot(str.alignErrInit)
    hold on
    plot(str.alignErr)
    legend('Init','Final')
    set(gca,'ylim',[0,5])
else
    plot(str.alignErr)
    legend('Final')
    set(gca,'ylim',[0,1])
end
xlabel('Spectral bins')
ylabel('Squared error (SO3)')
title('Alignment error')

subplot(Nrows,2,4)
plot(str.rc')
xlabel('Spectral bins')
ylabel('Rad')
title('Matching of spectral bins')


if Nrows == 4
    subplot(4,2,5)
    hold on

    if isfield(str,'ball')
        OA = str.ball.OA;
        mret = str.ball.mret;
        err = str.ball.err;
    else
        OA = str.OAball;
        mret = str.mret;
        err = str.err;
    end
    phi = unwrap(atan2(OA(2,:,:),OA(1,:,:)))/2;
    amp = sqrt(sum(OA.^2,1));
    plot(squeeze(amp.*cos(phi)),squeeze(amp.*sin(phi)),'-','LineWidth',2)

    phi = linspace(0,2*pi,201);
    amp = ones(1,201)*mret;
    plot(squeeze(amp.*cos(phi)),squeeze(amp.*sin(phi)),'-','LineWidth',2)

    axis equal
    xlabel('X')
    ylabel('Y')
    title(sprintf('Centered ball lens, err: %.3f',err))

    subplot(4,2,6)
    hold on
    phi = unwrap(atan2(OA(2,:),OA(1,:)));
    th = linspace(0,4*pi,numel(str.tissueRet));
    plot(th-phi - (mean(th-phi)),'LineWidth',2)
    xlabel('Eff rotation correction')
    ylabel('Rad')
    set(gca,'box','on','xlim',[1,numel(str.tissueRet)],'ylim',[-pi/2,pi/2])
    title('Vcorr eff')

    subplot(4,2,7)
    hold on
    plot(sqrt(sum(OA.^2,1)))
    plot([1,numel(str.tissueRet)],mret*[1,1])
    plot(str.tissueRet)
    xlabel('Rotation')
    ylabel('Ret')
    title('L2 ret')
    set(gca,'box','on','xlim',[1,numel(str.tissueRet)],'ylim',[0,4])

    subplot(4,2,8)
    hold on
    if numel(oaSheathStack)>0
      plot(linspace(.5,1,numel(oaSheathStack)).*cos(oaSheathStack),-linspace(.5,1,numel(oaSheathStack)).*sin(oaSheathStack),'.')
    end
    plot(str.tissueRet.*cos(str.tissueAngle),str.tissueRet.*sin(str.tissueAngle),'.-')
    plot((cos(str.sheathAngle(:))*[0,1])',(-sin(str.sheathAngle(:))*[0,1])','k')
    plot((cos(str.sheathAngle(:))*[-1,0])',(-sin(str.sheathAngle(:))*[-1,0])','--k')
    polar(linspace(0,2*pi,201),ones(1,201))%plot(ball.Vcomp)
    if ~isempty(oaxOffset)
        plot([0,cos(oaxOffset)],[0,-sin(oaxOffset)],'r')
    end
    axis equal
    xlabel('Q')
    ylabel('U')
    title(sprintf('Sheath,err %.3f',str.errSheathAngle))
    set(gca,'box','on')
end


