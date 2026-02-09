pathwayGenieApp.controller("pathwayGenieCtrl", ["PathwayGenieService", function(PathwayGenieService) {
	var self = this;
	
	self.toggleHelp = function() {
		return PathwayGenieService.toggleHelp();
	}
}]);