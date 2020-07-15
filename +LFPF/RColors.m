classdef RColors
    
    % A class to assign ROI colors in a consistent way
    
   properties
      CData
      ROIList %roi.Name
   end
   
    properties (Dependent)
        
    end
   
   
   %% ------------------
   methods
      function obj = RColors(ROIList)
          % receives a list of mouse V1/Thalamus ROIs and return proper
          % color for each. There is also subcolors for each for layer
          % colors
          if ~exist('ROIList','var') || isempty(ROIlist)
             ROIList = {'VISp','VISl','VISli','VISrl','VISal','VISpm','VISam','VISmma','LGd','LGv','LP'};
          end
          
          obj.ROIList = ROIList;
          for i = 1:numel(ROIList)
              % assign proper colors to ROIs
              switch ROIList{i}
                  case 'VISp'
                      Color = [0.43,0.25,0.63];
                  case 'VISl'
                      Color = [0.03,0.29,0.48];
                  case 'VISrl'
                      Color = [0.26,0.68,0.76];
                  case 'VISli'
                      Color = [0 .8 .8];
                  case 'VISal'
                      Color = [0.65,0.46,0.11];
                  case 'VISpm'
                      Color = [1 .7 .3];
                  case 'VISam'
                      Color = [0.8,0.11,0.11];
                  case 'VISmma'
                      Color = [0.40,0,0.12];
                  case 'LGv'
                      Color = [0.74,0.50,0.74];
                  case 'LGd'
                      Color = [0.97,0.51,0.75];
                  case 'LP'
                      Color = [0.3,0.67,0.29];
                  otherwise
                      Color = [.5 .5 .5];
                      
              end
              % define subcolors
              %Offset    = repmat(-.15:.05:.1,[3 1])';
              Offset    = repmat(-.25:.1:.25,[3 1])';
              SubColors = repmat(Color,[6 1])+Offset;
              SubColors = max(min(SubColors,1),0);
              %prepare output
              CData.(ROIList{i}).Color = Color; 
              CData.(ROIList{i}).SubColors = SubColors;
          end
          obj.CData = CData;
      end
      %%-----------------------------------------------
      function MColors = MatrixColors(obj,ROIlist,CType)
          % return colors in matrix/colormap format
          % ROIlist: the list of ROIs in the class, or a subset of it
          % CType, can be 'Color' or 'SubColors', for single color or 6 color per layer respectively
          
          if ~exist('ROIlist','var') || isempty(ROIlist)
              ROIlist = obj.ROIList;
          end
          
          if ~exist('CType','var') || isempty(CType)
              CType = 'Color';
          end
          
          if strcmpi(CType,'color')
              MColors = cellfun(@(x) obj.CData.(x).Color,ROIlist,'uni',false);
              MColors = cat(1,MColors{:});
          else
              MColors = cellfun(@(x) obj.CData.(x).SubColors,ROIlist,'uni',false);
              MColors = permute(cat(3,MColors{:}),[3 1 2]);
              
          end
          
          MColors = cat(1,MColors);
      end
   end
end