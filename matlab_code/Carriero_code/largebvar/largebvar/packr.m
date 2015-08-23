function datapackr=packr(DATAEVIEWS)
%remove NaNs
nanindex=isnan(DATAEVIEWS);

nanindex2=sum(nanindex,2)>0;

% size(nanindex2)
datapackr=DATAEVIEWS;
datapackr(nanindex2,:,:)=[];
