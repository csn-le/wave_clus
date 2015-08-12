function Plot_continuous_data(xf,handles,thr,thrmax)
% PLOT CONTINUOUS DATA
detect = handles.par.detection;
sr=handles.par.sr;
axes(handles.cont_data)
cla
hold on

plot((1:4:length(xf))/sr,xf(1:4:end))
%ylim([-thrmax thrmax])
switch detect
    case 'pos'
        line([0 length(xf)/sr], [thr thr],'color','r')
    case 'neg'
        line([0 length(xf)/sr], [-thr -thr],'color','r')
    case 'both'
        line([0 length(xf)/sr], [thr thr],'color','r')
        line([0 length(xf)/sr], [-thr -thr],'color','r')
end
axis([0 length(xf)/sr -thrmax/2 thrmax])
%xlabel('Time (sec)')
