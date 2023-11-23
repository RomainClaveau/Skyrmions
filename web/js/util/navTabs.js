var showNavTab = function(tab)
{
    document.querySelector("#tabs_nav span .selected").className = "";

    document.querySelector("#title_" + tab).className = "selected";

    document.querySelector(".tab.show").className = "tab hide";

    document.querySelector("#tab_" + tab).className = "tab show";
}

var showSectionTab = function(tab)
{
    document.querySelector("header span.selected").className = "";

    document.querySelector("#title_" + tab).className = "selected";

    document.querySelector(".sectionTab.show").className = "sectionTab hide";

    document.querySelector("#sectionTab_" + tab).className = "sectionTab show";

    if(tab == "3d")
    {
        document.querySelector("#chronology").style.display = "table";
        atomistic.constants.renderer.live = false;
        document.querySelector("section").style.height = "70%";
    }
    else
    {
        document.querySelector("#chronology").style.display = "none";
        atomistic.constants.renderer.live = false;
        document.querySelector("section").style.height = "80%";
        document.querySelector("#graph_Q_B_T").innerHTML = "Chargement du graphique...";
        setTimeout(function(){atomistic.studies.updateGraphics();}, 300);
    }

    document.querySelectorAll("#values_live input")[0].checked = atomistic.constants.renderer.live;
    
    atomistic.constants.renderer.currentTab = tab;
    
    if(tab == "table") atomistic.updateTable();
}