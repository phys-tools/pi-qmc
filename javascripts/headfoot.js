function writeHeader(relURL) {
    githubURL = "https://github.com/phys-tools/pi-qmc";
    document.write("<header>");
    document.write("  <h1>");
    document.write("    <a href='" + relURL + "index.html'>");
    document.write("      <img src='" + relURL + "images/pilogo.png'");
    document.write("           valign='bottom' height='48'/>");
    document.write("      &nbsp;");
    document.write("      pi-qmc");
    document.write("    </a>");
    document.write("  </h1>");
    document.write("  <p>Path integral quantum Monte Carlo</p>");
    document.write("  <p class='view'>");
    document.write("    <a href='" + githubURL + "'>");
    document.write("      View the Project on GitHub");
    document.write("      <small>pi-qmc</small>");
    document.write("    </a>");
    document.write("  </p>");
    document.write("  <ul>");
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
    document.write("</header>");
}

function writeFooter(relURL) {
    shumwayURL = "https://github.com/shumway";
    document.write("  <footer>");
    document.write("    <p>This project is maintained by");
    document.write("      <a href='" + shumwayURL + "'>shumway</a>");
    document.write("    </p>");
    document.write("    <p><small>");
    document.write("      Hosted by GitHub Pages &mdash;");
    document.write("      Theme by");
    document.write("      <a href='https://github.com/orderedlist'>");
    document.write("        orderedList</a>");
    document.write("    </small></p>"); 
    document.write("  </footer>");
}

