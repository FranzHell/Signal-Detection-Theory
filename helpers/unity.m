function unity(varargin);
% hn 08/23/06 
% accepts all the figure-property commands as 'plot.m' does
hold on;
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');

min_val = min([xlims,ylims]);
max_val = max([xlims,ylims]);

if nargin ==0
    plot([min_val max_val],[min_val max_val]);
else
%     plotstr = [];
%     for n = 1:length(varargin);
%         plotstr = [plotstr, varargin{n}];
%     end
%     plotstr
    plot([min_val max_val],[min_val max_val],varargin{:});
end
