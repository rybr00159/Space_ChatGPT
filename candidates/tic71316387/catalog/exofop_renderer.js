//
// Cell renderer - Errors
//
class ErrorCellRenderer {

  init(params) {
    // create the cell
    this.eGui = document.createElement('div');
    
    const data = params.api.getDisplayedRowAtIndex(params.rowIndex).data;
    const column = params.column;       // e.g. epoch
    const column_e = params.column_e;   // e.g. epoch_e
        
    this.eGui.innerHTML = `<span>`;
    
    if (data[params.column] || data[params.column] === 0) {
	    this.eGui.innerHTML += data[params.column];
		if (data[params.column_e]) {
			if (column == 'value') {
			this.eGui.innerHTML += " <span class=\"error3\">&plusmn; " + data[params.column_e] + "</span>";
			}
			else {
			this.eGui.innerHTML += " <div class=\"error3\">&plusmn; " + data[params.column_e] + "</div>";
			}
  		}
	}
	}
  getGui() {
    return this.eGui;
  }
}

//
// Cell renderer - Editing values
//
class EditCellRenderer {	
  init(params) {
    this.eGui = document.createElement('div');
    
    const data = params.api.getDisplayedRowAtIndex(params.rowIndex).data;    
    const edit_type = params.edit_type;   		// e.g. planet
    const edit_id = data[params.edit_id];       // e.g. 584
    const edit_user = data[params.edit_user];   // e.g. ciardi
    const my_user = params.my_user;   			// e.g. mcrane
    const my_status = params.my_status;   		// e.g. admin
//	const mag = params.mag;						// e.g. umag
	const mtabname = data[params.mtabname];		// e.g. umag
	const tid = params.tid;						// only for obs/files

    var url = '';
    this.eGui.innerHTML = `<span>`;

	//if ((edit_user == my_user || my_status == "admin") && (typeof edit_id !== 'undefined') ) {
	if ( (edit_user == my_user) && (typeof edit_id !== 'undefined') ) {

		if (edit_type == 'planet') url = 'edit_planet.php?pid=' + edit_id;
		else if (edit_type == 'planetparams') url = 'edit_planet_params.php?pid=' + edit_id;
		else if (edit_type == 'stellar') url = 'edit_stellar.php?sid=' + edit_id;
		else if (edit_type == 'scomp') url = 'edit_stellarcomp.php?scid=' + edit_id;
		else if (edit_type == 'mag') url = 'edit_mag.php?mag=' + mtabname + '&mid=' + edit_id;
		else if (edit_type == 'img') url = 'edit_img.php?imid=' + edit_id + '&tid=' + tid;
		else if (edit_type == 'spec') url = 'edit_spect.php?spid=' + edit_id + '&tid=' + tid;
		else if (edit_type == 'tser') url = 'edit_tseries.php?tsid=' + edit_id + '&tid=' + tid;
		else if (edit_type == 'file') url = 'edit_file.php?fid=' + edit_id + '&tid=' + tid;
				
		this.eGui.innerHTML += "<a class=\"over2\" href=\"#\" onClick=\"window.open('" + url + 
		"','_blank','width=800,height=1000,status=yes,resizable=yes,scrollbars=yes');return false;\">Edit</a>";
	}
    }
  getGui() {
    return this.eGui;
  }
}

//
// Cell renderer - View user info
//
class UserCellRenderer {
  init(params) {
  	this.eGui = document.createElement('div');
    
	const data = params.api.getDisplayedRowAtIndex(params.rowIndex).data;
	const user_id = data[params.user_id];
    const user_name = data[params.user_name];
    
    this.eGui.innerHTML = `<span>`;

    if (user_name == 'Exoplanet Archive' || user_name == 'TESS project') this.eGui.innerHTML += user_name;

    else {
	this.eGui.innerHTML += "<a class=\"over2\" href=\"#\"  onClick=\"window.open('view_obs.php?id=" +
			user_id + "', '_blank', 'width=500,height=500,status=yes,resizable=yes');return false;\">" +
			user_name + "</a>";
	}
	}
  getGui() {
    return this.eGui;
  }
}    

//
// Cell renderer - View tag info
//
class TagCellRenderer {
  init(params) {
  	this.eGui = document.createElement('div');
    
	const data = params.api.getDisplayedRowAtIndex(params.rowIndex).data;

    this.eGui.innerHTML = `<span>`;
	this.eGui.innerHTML += "<a class=\"over2\" href=\"#\"  onClick=\"window.open('view_tag.php?tag=" +
			params.value + "', '_blank', 'status=yes,resizable=yes');return false;\">" +
			params.value + "</a>";
	}
  getGui() {
    return this.eGui;
  }
}


//
// Cell renderer - Files
//
class FileCellRenderer {
  init(params) {
  	this.eGui = document.createElement('div');
    
	const data = params.api.getDisplayedRowAtIndex(params.rowIndex).data;
    const fid = data[params.fid];
    const fname = data[params.fname];
    const fext = data[params.fext];

    this.eGui.innerHTML = `<span>`;
    
    if (fext == 'jpg' || fext == 'jpeg' || fext == 'png' || fext == 'gif' || fext == 'pdf' || fext=='ps' ) {
    	this.eGui.innerHTML += "<a target=\"_blank\" href=\"get_file.php?id=" + fid + 
    	"\" data-container=\"body\" data-toggle=\"popover\" class=\"over2\" data-placement=\"right\" " +
    	"data-content=\"<img src=thumbtest.php?id=" + fid + "&width=225&height=225>\">" + fname + "</a>";
    }
    else {
    	this.eGui.innerHTML += "<a class=\"over2\" target=\"_blank\" href=\"get_file.php?id=" + fid + "\">" + fname + "</a>";
	}
	}
  getGui() {
    return this.eGui;
  }
}









