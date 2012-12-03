function writeHeader(relURL) {
    githubURL = "https://github.com/phys-tools/pi-qmc";
    document.write("<header>");
    writeTitle(relURL);
    writeViewProject(githubURL);
    writeGetCode(githubURL);
    writeNavigationLinks(relURL);
    document.write("</header>");
}

function writeTitle(relURL) {
    document.write("  <h1>");
    document.write("    <a href='" + relURL + "index.html'>");
    document.write("      <img src='" + relURL + "images/pilogo.png'");
    document.write("           valign='bottom' height='48'/>");
    document.write("      &nbsp;");
    document.write("      pi-qmc");
    document.write("    </a>");
    document.write("  </h1>");
    document.write("  <p>Path integral quantum Monte Carlo</p>");
}

function writeViewProject(githubURL) {
    document.write("  <p class='view'>");
    document.write("    <a href='" + githubURL + "'>");
    document.write("      View the Project on GitHub");
    document.write("      <small>pi-qmc</small>");
    document.write("    </a>");
    document.write("  </p>");
}

function writeGetCode(githubURL) {
    document.write("  <ul class='getCode'>");
    document.write("    <li>");
    document.write("      <a href='" + githubURL + "/zipball/master'>");
    document.write("        Download <strong>ZIP File</strong>");
    document.write("      </a>");
    document.write("    </li>");
    document.write("    <li>");
    document.write("      <a href='" + githubURL + "/tarball/master'>");
    document.write("        Download <strong>TAR Ball</strong>");
    document.write("      </a>");
    document.write("    </li>");
    document.write("    <li>");
    document.write("      <a href='" + githubURL + "'>");
    document.write("        View On <strong>GitHub</strong>");
    document.write("      </a>");
    document.write("    </li>");
    document.write("  </ul>");
}

function writeNavigationLinks(relURL) {
    shumwayURL = "http://shumway.physics.asu.edu";
    document.write("  <ul>");
    document.write("    <li>");
    document.write("      <a href='http://pi-qmc.readthedocs.org'>");
    document.write("        Documentation");
    document.write("      </a>");
    document.write("    </li>");
    document.write("    <li>");
    document.write("      <a href='" + relURL + "install/index.html'>");
    document.write("        Build / Install");
    document.write("      </a>");
    document.write("    </li>");
    document.write("    <li>");
    document.write("      <a href='" + relURL + "doxygen/index.html'>");
    document.write("        Doxygen");
    document.write("      </a>");
    document.write("    </li>");
    document.write("    <li>");
    document.write("      <a href='" + relURL + "team/people.html'>");
    document.write("        People");
    document.write("      </a>");
    document.write("    </li>");
    document.write("    <li>");
    document.write("      <a href='" + relURL + "team/groups.html'>");
    document.write("        Groups using pi-qmc");
    document.write("      </a>");
    document.write("    </li>");
    document.write("    <li>");
    document.write("      <a href='" + relURL + "team/funding.html'>");
    document.write("        Funding Sponsors");
    document.write("    <li>");
    document.write("      <a href='" + shumwayURL + "'>");
    document.write("        Shumway Research Group");
    document.write("      </a>");
    document.write("    </li>");
    document.write("  </ul>");
}

function writeFooter(relURL) {
    shumwayGitHubURL = "https://github.com/shumway";
    document.write("  <footer>");
    document.write("    <p>This project is maintained by");
    document.write("      <a href='" + shumwayGitHubURL + "'>shumway</a>");
    document.write("    </p>");
    document.write("    <p><small>");
    document.write("      Hosted by GitHub Pages &mdash;");
    document.write("      Theme by");
    document.write("      <a href='https://github.com/orderedlist'>");
    document.write("        orderedList</a>");
    document.write("    </small></p>"); 
    document.write("  </footer>");
}

