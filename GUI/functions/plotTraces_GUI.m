function fHandle = plotTraces_GUI(traces,spacing,vWindow,ax)
%fHandle = plotTraces_GUI(traces,spacing,vWindow,ax)
%   Produces a plot of each trace offset so that they aren't all
%   overlapping. Spacing sets how far apart the traces are (recommend
%   0.5-1) and vWindow sets the default vertical zoom as a multiple of
%   spacing. 
% if ~isempty(varargin) && strcmp(varargin{1}, 'plotSubWindow')
% fHandle = gcf;
% subplot(varargin{2}(1), varargin{2}(2), varargin{2}(3))
% 
% else
% fHandle = figure;
% 
% end

traces = bsxfun(@minus, traces, median(traces,2));
trace_visualization = bsxfun(@plus,traces,[0:size(traces,1)-1]'*spacing)';
plot(trace_visualization, 'Parent', ax);
ymin = min(trace_visualization(:));
ymax = max(trace_visualization(:));

% ylim([-spacing, vWindow*spacing]);
set(ax, 'XLim', [1, size(traces,2)])
set(ax, 'YLim', [ymin, ymax]);
set(ax, 'XTickMode', 'auto', 'YTickMode', 'auto');
axis(ax, 'normal');
% if ~isempty(varargin) && strcmp(varargin{1},'Label')
%     stack2Label = round(0.05*size(traces,2));
%     
%     for k=1:size(traces,1)
%         text(stack2Label,-spacing*(k-1)+.1*spacing,sprintf('%2.0f',varargin{2}(k)));
%     end
% end
