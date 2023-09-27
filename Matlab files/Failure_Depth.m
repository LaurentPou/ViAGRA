%close all;

if exist('frictionbis','var')
    
    if bool_plot_failure

        h=figure;
        plot((r_s(end)-r_s)/depth_unit,failure_radius(:,1));
        ylabel('Failure state (0: not reached, 1: reached)');
        xlabel(sprintf(['Depth (' depth_unit ')']));

    end
    
    min_failure = zeros(1,ff_size);
    max_failure = zeros(1,ff_size);
    min_cste = zeros(1,ff_size);
    max_cste = zeros(1,ff_size);
    min_failure2 = zeros(1,ff_size);
    max_failure2 = zeros(1,ff_size);
    min_cste2 = zeros(1,ff_size);
    max_cste2 = zeros(1,ff_size);
    
    for i = 1:ff_size
        % min_failure(i) = min(failure_radius_idx(failure_radius_idx(:,i) >
        % 0)); % bug incompréhensible...
        rr = 1;
        
        % Min surface failure determination
        while rr <= numel(r_s) && failure_radius_idx(rr,i) == 0
            rr = rr+1;
        end
        
        if rr < numel(r_s)
            min_failure(i) = failure_radius_idx(rr,i);
        else
            min_failure(i) = rr-1;
            min_cste(i) = -1;
        end
        
        if min_failure(i) == 0
            error('WTF');
        end
        
        % Max surface failure determination
        while rr <= numel(r_s)-1 && ...
                (failure_radius_idx(rr,i)+failure_radius_idx(rr+1,i) ~= 0)
            rr = rr+1;
        end
        
        if rr < numel(r_s)
            max_failure2(i) = failure_radius_idx(rr-1,i);
        else
            max_failure2(i) = rr;
            max_cste2(i) = 0;%-1;
        end
        
        % Min depth failure determination
        while rr <= numel(r_s)-1 && ...
                (failure_radius_idx(rr,i) + failure_radius_idx(rr+1,i) == 0)
            rr = rr+1;
        end
        
        if rr < numel(r_s)
            min_failure2(i) = failure_radius_idx(rr+1,i);
        else
            min_failure2(i) = rr;
            min_cste2(i) = 0;%-1;
        end
        
        
        % Max depth failure determination
        if ~isempty(failure_radius_idx(failure_radius_idx(:,i) > 0))
            max_failure(i) = max(failure_radius_idx(failure_radius_idx(:,i) > 0));
        else
            max_failure(i) = size(r_s);
            max_cste(i) = 0; %-1;
        end
        
        % Function continuity
        if max_failure2(i) == numel(r_s)
            max_failure2(i) = max_failure(i);
        end
        
        if min_failure2(i) == numel(r_s)
            min_failure2(i) = min_failure(i);
        end
        
    end
    
    if isequal(max_failure2,max_failure)
    
        h=figure;
        plot(frictionbis*180/pi,(r_s(end)-r_s(min_failure2))*depth_factor+min_cste2,'b');
        hold on;
        plot(frictionbis*180/pi,(r_s(end)-r_s(max_failure))*depth_factor+max_cste,'black');
        xlabel('Friction angle (°)');
        ylabel(sprintf(['Depth (' depth_unit ')']));
        legend('Max depth of failure','Min depth of failure');
        hold off;

        A_SavePlot(bool_save,h,sprintf(['Depth of failure ' model_name]));

        h=figure;
        semilogx((frictionbis*180/pi),(r_s(end)-r_s(min_failure2))*depth_factor+min_cste2,'b');
        hold on;
        semilogx((frictionbis*180/pi),(r_s(end)-r_s(max_failure))*depth_factor+max_cste,'black');
        xlabel('Friction angle (°)');
        ylabel(sprintf(['Depth (' depth_unit ')']));
        legend('Max depth of failure','Min depth of failure');
        hold off;
        % Update axis x labels
        New_XTickLabel = get(gca, 'xtick');
        set(gca,'XTickLabel',New_XTickLabel);
% 
%         % Updating for zoom
%         zoomH = zoom(gcf);
%         set(zoomH,'ActionPostCallback',{@A_Zoom_mypostcallback});

        A_SavePlot(bool_save,h,sprintf(['Depth of failure log x ' model_name]));

        h=figure;
        semilogy((frictionbis*180/pi),(r_s(end)-r_s(min_failure2))*depth_factor+min_cste2,'b');
        hold on;
        %semilogy((frictionbis*180/pi),(r_s(end)-r_s(max_failure))*depth_factor+max_cste,'black');
        xlabel('Friction angle (°)');
        ylabel(sprintf(['Depth (' depth_unit ')']));
        %legend('Max depth of failure','Min depth of failure');
        legend('Max depth of failure');
        hold off;
        % Update axis x labels
        New_YTickLabel2 = get(gca, 'ytick');
        set(gca,'YTickLabel',New_YTickLabel2);
% 
%         % Updating for zoom
%         zoomH = zoom(gcf);
%         set(zoomH,'ActionPostCallback',{@A_Zoom_mypostcallback});

        A_SavePlot(bool_save,h,sprintf(['Depth of failure log y ' model_name]));
    
    else
        
        h=figure;
        plot(frictionbis*180/pi,(r_s(end)-r_s(min_failure))*depth_factor+min_cste,'b');
        hold on;
        plot(frictionbis*180/pi,(r_s(end)-r_s(max_failure2))*depth_factor+max_cste2,'r');
        plot(frictionbis*180/pi,(r_s(end)-r_s(min_failure2))*depth_factor+min_cste2,'g');
        plot(frictionbis*180/pi,(r_s(end)-r_s(max_failure))*depth_factor+max_cste,'black');
        xlabel('Friction angle (°)');
        ylabel(sprintf(['Depth (' depth_unit ')']));
        legend('Max depth of deep failure','Min depth of deep failure', ...
            'Max depth of surface failure','Min depth of surface failure');
        hold off;

        A_SavePlot(bool_save,h,sprintf(['Depth of failure ' model_name]));

        h=figure;
        semilogx((frictionbis*180/pi),(r_s(end)-r_s(min_failure))*depth_factor+min_cste,'b');
        hold on;
        semilogx((frictionbis*180/pi),(r_s(end)-r_s(max_failure2))*depth_factor+max_cste2,'r');
        semilogx((frictionbis*180/pi),(r_s(end)-r_s(min_failure2))*depth_factor+min_cste2,'g');
        semilogx((frictionbis*180/pi),(r_s(end)-r_s(max_failure))*depth_factor+max_cste,'black');
        xlabel('Friction angle (°)');
        ylabel(sprintf(['Depth (' depth_unit ')']));
        legend('Max depth of deep failure','Min depth of deep failure', ...
            'Max depth of surface failure','Min depth of surface failure');
        hold off;
        % Update axis x labels
        New_XTickLabel = get(gca, 'xtick');
        set(gca,'XTickLabel',New_XTickLabel);
% 
%         % Updating for zoom
%         zoomH = zoom(gcf);
%         set(zoomH,'ActionPostCallback',{@A_Zoom_mypostcallback});

        A_SavePlot(bool_save,h,sprintf(['Depth of failure log x ' model_name]));

        h=figure;
        semilogy((frictionbis*180/pi),(r_s(end)-r_s(min_failure))*depth_factor+min_cste,'b');
        hold on;
        semilogy((frictionbis*180/pi),(r_s(end)-r_s(max_failure2))*depth_factor+max_cste2,'r');
        semilogy((frictionbis*180/pi),(r_s(end)-r_s(min_failure2))*depth_factor+min_cste2,'g');
        %semilogy((frictionbis*180/pi),(r_s(end)-r_s(max_failure))*depth_factor+max_cste,'black');
        xlabel('Friction angle (°)');
        ylabel(sprintf(['Depth (' depth_unit ')']));
        %legend('Max depth of deep failure','Min depth of deep failure', ...
        %    'Max depth of surface failure','Min depth of surface failure');
        legend('Max depth of deep failure','Min depth of deep failure', ...
            'Max depth of surface failure');
        hold off;
        % Update axis x labels
        New_YTickLabel2 = get(gca, 'ytick');
        set(gca,'YTickLabel',New_YTickLabel2);
% 
%         % Updating for zoom
%         zoomH = zoom(gcf);
%         set(zoomH,'ActionPostCallback',{@A_Zoom_mypostcallback});

        A_SavePlot(bool_save,h,sprintf(['Depth of failure log y ' model_name]));

    end
        
else
    
    if bool_plot_failure

        h=figure;
        plot((r_s(end)-r_s)*depth_factor,failure_radius_rel);
        ylabel('Failure state (0: not reached, 1: reached)');
        xlabel(sprintf(['Depth (' depth_unit ')']));
    
        h=figure;
        plot(failure_radius_rel,r_s*depth_factor);
        xlabel('Failure state (0: not reached, 1: reached)');
        ylabel(sprintf(['Radius (' depth_unit ')']));

    end
    
end