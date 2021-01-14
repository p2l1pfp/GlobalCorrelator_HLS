library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity calo_router is
    port(
            ap_clk : IN STD_LOGIC;
            enabled : IN STD_LOGIC;
            newevent : IN STD_LOGIC;
            links_in      : IN particles(NCALOSECTORS*NCALOFIBERS-1 downto 0);
            fifo_in       : OUT particles(NCALOSECTORS*NCALOFIFOS-1 downto 0);
            fifo_in_write : OUT std_logic_vector(NCALOSECTORS*NCALOFIFOS-1 downto 0);
            fifo_in_roll  : OUT std_logic_vector(NCALOSECTORS*NCALOFIFOS-1 downto 0)
    );
end calo_router;

architecture Behavioral of calo_router is
begin

    link2fifo : process(ap_clk)
        variable inext, ilink, ififo : natural;
    begin
        if rising_edge(ap_clk) then
            for isec in 0 to NCALOSECTORS-1 loop
                -- derive indices of previous and next sector
                if isec <= 1 then 
                    inext := isec+1; 
                else 
                    inext := 0; 
                end if;
                -- first region, that has only one calo sector contributing, if |phi| <= PHI_HALFWIDTH
                for ifib in 0 to NCALOFIBERS-1 loop
                    ilink := isec*NCALOFIBERS+ifib; ififo := isec*NCALOFIFOS+ifib;
                    fifo_in(ififo) <= links_in(ilink);
                    fifo_in_roll(ififo) <= newevent;
                    if enabled = '0' or links_in(ilink).pt = 0 then
                        fifo_in_write(isec*NCALOFIFOS+ifib) <= '0';
                    else
                        if links_in(ilink).phi <= PHI_HALFWIDTH_POS and
                           links_in(ilink).phi >= PHI_HALFWIDTH_NEG then
                            fifo_in_write(ififo) <= '1';
                        else
                            fifo_in_write(ififo) <= '0';
                        end if;
                    end if;
                end loop; 
                -- second region 
                ---- take from this sector if phi >= PHI_MARGIN_POS, shift back by 1 region size
                for ifib in 0 to NCALOFIBERS-1 loop
                    ilink := isec*NCALOFIBERS+ifib; ififo := isec*NCALOFIFOS+NCALOFIBERS+ifib;
                    fifo_in(ififo).pt   <= links_in(ilink).pt;
                    fifo_in(ififo).eta  <= links_in(ilink).eta;
                    fifo_in(ififo).phi  <= links_in(ilink).phi - PHI_SHIFT;
                    fifo_in(ififo).rest <= links_in(ilink).rest;
                    fifo_in_roll(ififo) <= newevent;
                    if enabled = '0' or links_in(ilink).pt = 0 then
                        fifo_in_write(ififo) <= '0';
                    else
                        if links_in(ilink).phi >= PHI_MARGIN_POS then
                            fifo_in_write(ififo) <= '1';
                        else
                            fifo_in_write(ififo) <= '0';
                        end if;
                    end if;
                end loop; 
                ---- take from the next sector if phi <= PHI_CALOEDGE_NEG, shift up by 1 sector, down by 1 region
                for ifib in 0 to NCALOFIBERS-1 loop
                    ilink := inext*NCALOFIBERS+ifib; ififo := isec*NCALOFIFOS+2*NCALOFIBERS+ifib;
                    fifo_in(ififo).pt   <= links_in(ilink).pt;
                    fifo_in(ififo).eta  <= links_in(ilink).eta;
                    fifo_in(ififo).phi  <= links_in(ilink).phi + PHI_CALOSHIFT1;
                    fifo_in(ififo).rest <= links_in(ilink).rest;
                    fifo_in_roll(ififo) <= newevent;
                    if enabled = '0' or links_in(ilink).pt = 0 then
                        fifo_in_write(ififo) <= '0';
                    else
                        if links_in(ilink).phi <= PHI_CALOEDGE_NEG then
                            fifo_in_write(ififo) <= '1';
                        else
                            fifo_in_write(ififo) <= '0';
                        end if;
                    end if;
                end loop;
                -- third region 
                ---- take from this sector if phi >= PHI_CALOEDGE_POS
                for ifib in 0 to NCALOFIBERS-1 loop
                    ilink := isec*NCALOFIBERS+ifib; ififo := isec*NCALOFIFOS+3*NCALOFIBERS+ifib;
                    fifo_in(ififo).pt   <= links_in(ilink).pt;
                    fifo_in(ififo).eta  <= links_in(ilink).eta;
                    fifo_in(ififo).phi  <= links_in(ilink).phi - PHI_CALOSHIFT1;
                    fifo_in(ififo).rest <= links_in(ilink).rest;
                    fifo_in_roll(ififo) <= newevent;
                    if enabled = '0' or links_in(ilink).pt = 0 then
                        fifo_in_write(ififo) <= '0';
                    else
                        if links_in(ilink).phi >= PHI_CALOEDGE_POS then
                            fifo_in_write(ififo) <= '1';
                        else
                            fifo_in_write(ififo) <= '0';
                        end if;
                    end if;
                end loop; 
                ---- take from the next sector if phi <= PHI_MARGIN_NEG
                for ifib in 0 to NCALOFIBERS-1 loop
                    ilink := inext*NCALOFIBERS+ifib; ififo := isec*NCALOFIFOS+4*NCALOFIBERS+ifib;
                    fifo_in(ififo).pt   <= links_in(ilink).pt;
                    fifo_in(ififo).eta  <= links_in(ilink).eta;
                    fifo_in(ififo).phi  <= links_in(ilink).phi + PHI_SHIFT;
                    fifo_in(ififo).rest <= links_in(ilink).rest;
                    fifo_in_roll(ififo) <= newevent;
                    if enabled = '0' or links_in(ilink).pt = 0 then
                        fifo_in_write(ififo) <= '0';
                    else
                        if links_in(ilink).phi <= PHI_MARGIN_NEG then
                            fifo_in_write(ififo) <= '1';
                        else
                            fifo_in_write(ififo) <= '0';
                        end if;
                    end if;
                end loop;
            end loop;
        end if;
    end process link2fifo;


end Behavioral;
