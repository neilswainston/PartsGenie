resultApp.factory("ResultService", ["$http", "$rootScope", "$window", "ErrorService", "ProgressService", function($http, $rootScope, $window, ErrorService, ProgressService) {
	var obj = {};
	obj.results = null;
	obj.response = {"update": {}};
	resultsSaved = false;
	
	var jobIds = [];
	var jobId = null;
	
	obj.setResults = function(res) {
		resultsSaved = false;
		obj.results = res;
	};
	
	obj.appendResults = function(res) {
		resultsSaved = false;
		
		if(!obj.results) {
			obj.results = [];
		}
		
		// Note: res has changed from array to single object:
		obj.results.push(res);
	};
	
	obj.resultsSaved = function() {
		return resultsSaved;
	}
	
	obj.cancel = function() {
		return PathwayGenieService.cancel(jobId);
	};
	
	obj.update = function() {
		return obj.response.update;
	};
	
	obj.listen = function() {
		if(jobIds.length == 0) {
			resultsSaved = true;
			return;
		}
		
		jobId = jobIds[0];
		var source = new EventSource("/progress/" + jobId);

		source.onmessage = function(event) {
			obj.response = JSON.parse(event.data);
			status = obj.response.update.status;
			
			if(status == "cancelled" || status == "error" || status == "finished") {
				source.close();

				if(status == "finished") {
					for (i = 0; i < obj.response.result.length; i++) {
						ice_ids = obj.response.result[i];
						
						for(var key in ice_ids) {
							link = ice_ids[key].link;
							
							if(obj.results[i].links.indexOf(link) == -1 ) {
								obj.results[i].links.push(link);
							}
						}
						obj.results[i].ice_ids = ice_ids;
					}
				}
				
				jobIds.splice(0, 1);
				obj.listen();
			}
			
			$rootScope.$apply();
		};
		
		source.onerror = function(event) {
			source.close();
			jobIds.splice(0, 1);
			listen();
			onerror(event.message);
		}
	};

	obj.exportOrder = function() {
		$http.post("/export", {"designs": obj.results}).then(
				function(resp) {
					var newWindow = $window.open();
					newWindow.location = resp.data.path;
				},
				function(errResp) {
					ErrorService.open(errResp.data.message);
				});
	};

	return obj;
}]);