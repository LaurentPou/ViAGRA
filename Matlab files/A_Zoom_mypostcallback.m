function A_Zoom_mypostcallback(~,~)

% Updates ticks after zooming on plot (useful for semilog plots)

New_XTickLabel = get(gca,'xtick');
set(gca,'XTickLabel',New_XTickLabel);

end

% Does not work when translating !

% Example of use with semilog plots
% figure;
% semilogx((frictionbis*180/pi),r_s(end)-r_s(min_failure)+min_cste,'b');
% hold on;
% semilogx((frictionbis*180/pi),r_s(end)-r_s(max_failure2)+max_cste2,'r');
% semilogx((frictionbis*180/pi),r_s(end)-r_s(min_failure2)+min_cste2,'g');
% semilogx((frictionbis*180/pi),r_s(end)-r_s(max_failure)+max_cste,'black');
% xlabel('Friction angle (°)');
% ylabel('Depth (m)');
% legend('Max depth of deep failure','Min depth of deep failure', ...
%     'Max depth of surface failure','Min depth of surface failure');
% hold off;
% % Update axis x labels
% New_XTickLabel = get(gca, 'xtick');
% set(gca,'XTickLabel',New_XTickLabel);
% 
% % Updating for zoom
% zoomH = zoom(gcf);
% set(zoomH,'ActionPostCallback',{@A_Zoom_mypostcallback});