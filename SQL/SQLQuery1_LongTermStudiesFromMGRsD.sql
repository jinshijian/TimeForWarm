
WITH CET (StudyNumber, MeasureYear) AS 
(
	SELECT DISTINCT StudyNumber, Measure_year FROM [dbo].[SRDBMonthly]
	GROUP BY [StudyNumber], Measure_year	
)

SELECT StudyNumber, Count(StudyNumber) AS obs_yr FROM CET GROUP BY StudyNumber
ORDER BY obs_yr DESC